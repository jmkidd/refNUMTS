import sys
import genutils3
import os
import numpy as np
import shutil
import argparse

#### START MAIN PROGRAM ########
parser = argparse.ArgumentParser(description='find numts in assembly')

parser.add_argument('--sourcename', type=str,help='name of source genome',required=True)
parser.add_argument('--sourcefasta', type=str,help='source genome fasta',required=True)
parser.add_argument('--sourceregions', type=str,help='NUMT ass_table.txt file for source',required=True)

parser.add_argument('--targetname', type=str,help='name of target genome',required=True)
parser.add_argument('--targetfasta', type=str,help='target genome fasta',required=True)
parser.add_argument('--outdirbase', type=str,help='output directory base dir',required=True)


parser.add_argument('--mitoref', type=str,help='mitochondria fasta with .fai',required=True)
parser.add_argument('--mitorefdouble', type=str,help='mitochondria fasta concatenated 2x with .fai',required=True)

args = parser.parse_args()

sourceGenomeName = args.sourcename
sourceGenomeFasta = args.sourcefasta

targetGenomeName = args.targetname
targetGenomeFasta = args.targetfasta
outDirBase = args.outdirbase


# get fasta of each assembled numt in source
sourceRegions = args.sourceregions

# setup some needed info and default values
# value of intervals to add on each end for searching
intPad = 1000
minFlankCovBp = 200
minIntCovFrac = 0.8

cToDo = []
for i in range(1,39):
	c = 'chr' + str(i)
	cToDo.append(c)
cToDo.append('chrX')


if os.path.isdir(outDirBase) is False:
    print('ERROR! %s not found' % outDirBase)
    sys.exit()
if outDirBase[-1] != '/':
    outDirBase += '/'

outDir = outDirBase + sourceGenomeName + '_in_' + targetGenomeName
cmd = 'mkdir %s ' % outDir
if os.path.isdir(outDir) is False:
    genutils3.runCMD(cmd)
outDir += '/'


print('reading in assembled numt regions from %s' % sourceRegions)

assemIntsToDo = {}
inFile = open(sourceRegions,'r')
for line in inFile:
    if line[0] == '#':
        continue
    line = line.rstrip()
    line = line.split()
    assemID = line[0]
    assemInt = line[1]
    assemLen = int(line[13])
    hspID = line[2]
    chrm = line[3]
    if chrm not in cToDo:
        continue
    assemIntsToDo[assemID] = [assemInt,assemLen]
inFile.close()

print('have %i assemInts to do' % len(assemIntsToDo))


assemSeqsFA = outDir + 'assemSeqs_pad.fa'
tmpFA = outDir + 'tmp.fa'
outFile = open(assemSeqsFA,'w')

alignDecision = {}
for assemID in assemIntsToDo.keys():
    alignDecision[assemID] = []

#    print(assemID,assemIntsToDo[assemID])
    p = assemIntsToDo[assemID][0].split(':')
    cName = p[0]
    p = p[1].split('-')
    start = int(p[0])
    end = int(p[1])    
#    print(cName,start,end)
    newStart = start - intPad
    newEnd = end + intPad
    
    extReg = '%s:%i-%i' % (cName,newStart,newEnd)
    
    extCmd = 'samtools faidx %s %s > %s' % (sourceGenomeFasta,extReg,tmpFA)
#    print(extCmd)
    genutils3.runCMD(extCmd)
    
    extSeq = genutils3.read_fasta_file_to_dict(tmpFA)
    seq = extSeq[extReg]['seq']
    seq = genutils3.add_breaks_to_line(seq,n=100)
    outFile.write('>%s\n%s\n' % (assemID+'_ext',seq))
    
outFile.close()

# now, search against target genome using minimap2
# should be minimap2/2.26

outPAFFile = outDir + sourceGenomeName + '_vs_' + targetGenomeName + '.paf'
cmd = 'minimap2 -x asm20 -c -o %s %s %s' % (outPAFFile,targetGenomeFasta,assemSeqsFA)
print(cmd)
genutils3.runCMD(cmd)

# now go through each and parse out potential hits...
# require their to be a single hit the spans across??

inFile = open(outPAFFile,'r')
for line in inFile:
    line = line.rstrip()
    line = line.split()
    pafFields = genutils3.parse_paf_line(line)
    cigarStr = line[-1]
    cigarStr = cigarStr.split(':')
    if cigarStr[0] != 'cg':
        print('error! CIGAR not found!')
        print(line)
        sys.exit()
    cigarStr = cigarStr[2]
    cigarExpand = genutils3.expand_cigar(cigarStr)
    
    qNameOrg = pafFields['qName'].replace('_ext','')
    # setup align mask
    qCovg = np.zeros(pafFields['qLen']+1,dtype=int)
    
    if pafFields['tName'] == 'chrM':
        continue
    
    tHit = pafFields['tName'] + ':' + str(pafFields['tStart'])+ '-' + str(pafFields['tEnd'])
        
    
    currentPos = 0
    if pafFields['strand'] == '+':
        currentPos = pafFields['qStart'] + 1 # do it 1 based
        for ops in cigarExpand:
            if ops[1] == 'M':
                for i in range(ops[0]):
                    qCovg[currentPos] = 1
                    currentPos += 1
            elif ops[1] == 'I':
                for i in range(ops[0]):
                    currentPos += 1
            elif ops[1] == 'D':
                continue
            else:
                print('what?',ops)
                sys.exit()
            
        
    elif pafFields['strand'] == '-':
        currentPos = pafFields['qEnd']
        for ops in cigarExpand:
            if ops[1] == 'M':
                for i in range(ops[0]):
                    qCovg[currentPos] = 1
                    currentPos -= 1
            elif ops[1] == 'I':
                for i in range(ops[0]):
                    currentPos -= 1
            elif ops[1] == 'D':
                continue
            else:
                print('what?',ops)
                sys.exit()
        

    leftFlankStart = 1
    leftFlankEnd = leftFlankStart+intPad-1
        
    rightFlankEnd = pafFields['qLen']
    rightFlankStart = rightFlankEnd - intPad + 1
        
    midStart = leftFlankEnd + 1
    midEnd = rightFlankStart - 1
        
    covLeft = qCovg[leftFlankStart:leftFlankEnd+1].sum()
    covMid = qCovg[midStart:midEnd+1].sum()
    covRight = qCovg[rightFlankStart:rightFlankEnd+1].sum()
    
    # require at least 200 bp of left and right and at least 80% of midd
    midLen = midEnd-midStart + 1
    midPer = covMid/midLen
    
    covConclusion = '?'
    
    if covLeft >= minFlankCovBp and covRight >= minFlankCovBp and midPer >= minIntCovFrac:
        covConclusion = 'PRESENT'
    elif covLeft >= minFlankCovBp and covRight >= minFlankCovBp:
        covConclusion = 'FLANKS_ONLY'
    elif covLeft >= minFlankCovBp :
        covConclusion = 'LEFT_ONLY'
    elif covRight >= minFlankCovBp :
        covConclusion = 'RIGHT_ONLY'
    else:
        covConclusion = 'NOT_FOUND'
        
    
    if covConclusion != 'NOT_FOUND':
        alignDecision[qNameOrg].append([covConclusion,tHit])        
inFile.close()

outMapSummaryFile = outPAFFile + '.summary.txt'
outFile = open(outMapSummaryFile,'w')

nl = ['assemID','assemInterval','assemLen','searchResults']
nl = '\t'.join(nl) + '\n'
outFile.write(nl)

for assemID in assemIntsToDo.keys():
    nl = [assemID,assemIntsToDo[assemID][0],str(assemIntsToDo[assemID][1])]
    
    hits = []
    if len(alignDecision[assemID]) == 0:
        hits = ['NO_HITS']
    else:
        for i in alignDecision[assemID]:
            j = ','.join(i)
            hits.append(j)
    hits = ';'.join(hits)
    nl.append(hits)
    nl = '\t'.join(nl) + '\n'
    outFile.write(nl)    
outFile.close()

print('hits summary written to',outMapSummaryFile)

