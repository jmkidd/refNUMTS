import sys
import os
import subprocess
import shutil
import argparse


#####################################################################
# Helper function to run commands, handle return values and print to log file
def runCMD(cmd):
    val = subprocess.Popen(cmd, shell=True).wait()
    if val == 0:
        pass
    else:
        print('command failed')
        print(cmd)
        sys.exit(1)
#############################################################################        
# setup paths to default programs to use and checks for required programs
def check_prog_paths(myData):        
    print('\nChecking for required programs...\n')
    
    for p in ['samtools','sortBed','blastn']:
        if shutil.which(p) is None:
            s = p + ' not found in path! please fix (module load?)'
            print(s, flush=True)
            sys.exit()
        else:
            print(p,shutil.which(p),flush=True)
            
#####################################################################
def run_get_genomes_in_assem(myData):
    # This is the main driver script that runs all of the analysis steps
    check_prog_paths(myData)
    print(myData['genomeName'], myData['genomeFA'])
    genomeOutDir = myData['genomeOutDirBase'] + myData['genomeName']
    if os.path.isdir(genomeOutDir) is False:
        print('make genome dir! %s' % genomeOutDir)
        cmd = 'mkdir %s ' % genomeOutDir
        print(cmd)
        runCMD(cmd)
    myData['genomeOutDir'] = genomeOutDir + '/'
    run_mito_blast2seq(myData)
    parse_blast_output(myData)    
    merge_blast_output(myData)
    align_numt_fragments(myData)
    make_combined_numt_fragments_table(myData)
#####################################################################
def run_mito_blast2seq(myData):
    # this runs blast2seq with blastn of the mito fasta
    # against the genome fasta
    print('in run_mito_blast2seq')
    myData['mitoBLASTout'] = myData['genomeOutDir'] + myData['genomeName'] + '.blast.out'
    if os.path.isfile(myData['mitoBLASTout']) is False:
        print('running blast!')
        
        # blastn options and alignment criteria from Simone et al. https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-12-517
        cmd = 'blastn -task blastn -reward 2 -penalty -3 -gapopen 5 -gapextend 2 -evalue 0.001 -outfmt 7 -dust no -query %s -subject %s -out %s \n' % (myData['mitoRefFA'], myData['genomeFA'],myData['mitoBLASTout'])
        print(cmd)
        runCMD(cmd)
#####################################################################
def parse_blast_output(myData):
    # performs initial parse of blast2seq output
    myData['blastParse'] = myData['mitoBLASTout'] + '.parse'
    myData['blastParseSort'] = myData['blastParse'] + '.sort'
    
    if os.path.isfile(myData['blastParseSort']) is True:
        print('skipping, already parsed and sorted!')
        return
    
    inFile = open(myData['mitoBLASTout'],'r')
    outFile = open(myData['blastParse'],'w')
    for line in inFile:
        if line[0] == '#':
            continue
        line = line.rstrip()
        line = line.split()
        qName = line[0]
        sName = line[1]
        pID = line[2]
        qStart = line[6]
        qEnd = line[7]
        sStart = line[8]
        sEnd = line[9]
        eVal = line[10]
        
        # filter        
        if float(eVal) > myData['minEval']: # skip min evals
            continue
        if sName in myData['chromsToSkipHits']:
            continue

        sDir = '?'
        if int(sStart) < int(sEnd):
            sDir = '+'
        else:
            sDir = '-'
            t = sStart
            sStart = sEnd
            sEnd = t        
        nl = [sName,sStart,sEnd,qName,qStart,qEnd,pID,eVal,sDir]
        
        nl = '\t'.join(nl) + '\n'
        outFile.write(nl)

    inFile.close()
    outFile.close()
    
    # sort
    cmd = 'sortBed -i %s > %s ' % (myData['blastParse'],myData['blastParseSort'])
    print(cmd)
    runCMD(cmd)
#####################################################################
def merge_blast_output(myData):
    # merges initial hits into assembled loci
    # based on distance on the genome and on the mito chromosome
    print('starting merge!')
    print('merge delta threshold is',myData['mergeDelta'])
    myData['blastParseSortMerged'] = myData['blastParseSort'] + '.merged'
    if os.path.isfile(myData['blastParseSortMerged']) is True:
        print('skipping merge, already done!')
        return
        
    hitsPerChrom = {}
    inFile = open(myData['blastParseSort'],'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        t = line[0]
        tS = int(line[1])
        tE = int(line[2])
        q = line[3]
        qS = int(line[4])
        qE = int(line[5])
        orient = line[8]
       
        nl = [t,tS,tE,q,qS,qE,orient]
        if t not in hitsPerChrom:
            hitsPerChrom[t] = []
        hitsPerChrom[t].append(nl)
    inFile.close()
    print(len(hitsPerChrom))
    for tName in hitsPerChrom:
        print('starting %s %i' % (tName,len(hitsPerChrom[tName])))
        
        moreToDo = True
        while moreToDo is True:
            if len(hitsPerChrom[tName]) == 1: # nothing to merge
                moreToDo = False
                break
            numMerges = 0
            for i in range(0,len(hitsPerChrom[tName])-1):
                if numMerges > 0: # as list index might be off..
                    break
                for j in range(i+1,len(hitsPerChrom[tName])):
                    # check to see if can do a merge
                    # det delta along t and q
                    tDelt = hitsPerChrom[tName][j][1] - hitsPerChrom[tName][i][2]

                    # query coords might not be in order                    
                    qHits = [ [hitsPerChrom[tName][j][4],hitsPerChrom[tName][j][5]],
                               [hitsPerChrom[tName][i][4],hitsPerChrom[tName][i][5]] ]
                    qHits.sort()
                    qDelt = qHits[1][0] - qHits[0][1]
                    
                    # check to see if we have a hit that is near 
                    # within 2000 if go around the circle...
                    diffFromStart = qHits[0][0] + (myData['mitoRefLen'] - qHits[1][1])
                    if diffFromStart < myData['mergeDelta']:
                        qDelt = diffFromStart
                    
                    # diff orients
                    if hitsPerChrom[tName][i][6] != hitsPerChrom[tName][j][6]:
                       continue
                    
                    
                    if tDelt > myData['mergeDelta']:
                       continue
                    
                    if tDelt <= myData['mergeDelta'] and qDelt <= myData['mergeDelta']:
                        numMerges += 1
                        # do merge                        
                        hitsPerChrom[tName][i][2] = max(hitsPerChrom[tName][j][2],hitsPerChrom[tName][i][2])
                        hitsPerChrom[tName][i][4] = qHits[0][0]                        
                        hitsPerChrom[tName][i][5] = qHits[1][1]                        
                        
                        hitsPerChrom[tName].pop(j)
                        break
                if numMerges > 0:
                    moreToDo = True
                else:
                    moreToDo = False
        
        print('ending %s %i' % (tName,len(hitsPerChrom[tName])))
        outFile = open(myData['blastParseSortMerged'],'w')
        for tName in hitsPerChrom:
            for i in hitsPerChrom[tName]:
                nl = [str(j) for j in i]
                nl = '\t'.join(nl) + '\n'
                outFile.write(nl)        
        outFile.close()

#####################################################################
def align_numt_fragments(myData):
    # align the numt fragmnets using blast2seq again to get new coordiates
    # and to get the sequence identity
    # this searches againts a fasta that has the mito sequence repeated
    # so that hsps that cross the circular boundary can be found
    # this allows for a better HSP count and for estimate of sequence identity
    
    print('running align_numt_fragments...')
    myData['fragAlignDir'] = myData['genomeOutDir'] + 'fragalign'
    if os.path.isdir(myData['fragAlignDir'] ) is False:
        print('make frag align dir! %s' % myData['fragAlignDir'] )
        cmd = 'mkdir %s ' % myData['fragAlignDir'] 
        runCMD(cmd)
    myData['fragAlignDir'] += '/'
    
    inFile = open(myData['blastParseSortMerged'],'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        
        extCoords = line[0] + ':' + line[1] + '-' + line[2]
        extFileName = myData['fragAlignDir'] + extCoords + '.fa'
        extBlast = extFileName + '.out'
        extBlastParse = extBlast + '.parse'

        
        cName = line[0]
        cExtStart = int(line[1])
        
        segLen = int(line[2]) - int(line[1]) + 1
        
        if os.path.isfile(extFileName) is True and os.path.isfile(extBlast) is True and os.path.isfile(extBlastParse) is True:
            continue
        
        # extract fragments and search with blast2seq                
        extCmd = 'samtools faidx %s %s > %s ' % (myData['genomeFA'],extCoords,extFileName)
        runCMD(extCmd)
                
        cmd = 'blastn -task blastn -reward 2 -penalty -3 -gapopen 5 -gapextend 2 -evalue 0.001 -outfmt 7 -dust no -query %s -subject %s -out %s \n' % (myData['mitoRefFADouble'], extFileName,extBlast)
        runCMD(cmd)

        
        hits = read_in_blast_hits(extBlast)
        # go through and update them all to have proper chrom coords based on extraction
        for i in range(len(hits)):
            hits[i][0] = cName
            hits[i][1] = hits[i][1] + cExtStart -1
            hits[i][2] = hits[i][2] + cExtStart -1        
        
        # remove hits that are entirely in 2nd copy of mito
        toRemove = []
        for i in range(len(hits)):
            if hits[i][4] >  myData['mitoRefLen'] and hits[i][5] > myData['mitoRefLen']:
                toRemove.append(i)
        for i in sorted(toRemove,reverse=True):
            del hits[i]
        
        
        # now, remove ones that fully overlap, don't want these sub hits
        moreToDo = True
        while moreToDo is True:
            if len(hits) == 1: # nothing to merge
                moreToDo = False
                break
            numMerges = 0
            for i in range(0,len(hits)-1):
                if numMerges > 0: # as list index might be off..
                    break
                for j in range(i+1,len(hits)):
                    if hits[j][1] >= hits[i][1] and hits[j][2] <= hits[i][2]:
                       numMerges +=1
                       hits.pop(j)
                       break
            if numMerges == 0:
                moreToDo = False
        
        
       # fix names and coordinates
        for i in range(len(hits)):       
            hits[i][3] = myData['mitoFAName'] 
            if hits[i][4] < myData['mitoRefLen'] and hits[i][5] <= myData['mitoRefLen']: # nothing to do
                continue
                
            start1 = hits[i][4]
            end1 = myData['mitoRefLen']
            start2 = 1
            end2 = hits[i][5] - myData['mitoRefLen'] 
            
            # convert to starts and ends blocks when cross boundary
            hits[i][4] = '%i,%i' % (start1,start2)
            hits[i][5] = '%i,%i' % (end1,end2)            

        # sort hits
        hits.sort(key=lambda x: x[1])

        outFile = open(extBlastParse,'w')
        for i in hits:
            nl = [str(j) for j in i]
            nl = '\t'.join(nl) + '\n'
            outFile.write(nl)
        outFile.close()     
    inFile.close() # close of all merged/assembled loci
#####################################################################
def read_in_blast_hits(fileName):
    # read in blast hits into list of lists
    hits = []
    inFile = open(fileName,'r')
    for line in inFile:
        if line[0] == '#':
            continue
        line = line.rstrip()
        line = line.split()
        qName = line[0]
        sName = line[1]
        pID = float(line[2])
        qStart = int(line[6])
        qEnd = int(line[7])
        sStart = int(line[8])
        sEnd = int(line[9])
        eVal = float(line[10])
        sDir = '?'
        if sStart < sEnd:
            sDir = '+'
        else:
            sDir = '-'
            t = sStart
            sStart = sEnd
            sEnd = t        
        nl = [sName,sStart,sEnd,qName,qStart,qEnd,pID,eVal,sDir]
        hits.append(nl)
    return hits
#####################################################################
def make_combined_numt_fragments_table(myData):
    # makes summary table based on alignment to the doubled mito fasta
    # adds assembly and HSP ids to the table
    print('running make_combined_numt_fragments_table...')    
    myData['assemTable'] = myData['blastParseSortMerged'] + '.assem_table.txt'
    inFile = open(myData['blastParseSortMerged'],'r')
    outFile = open(myData['assemTable'],'w')
    
    nl = ['#assemID','assemInterval','hspID','genomeChr','genomeStart','genomeEnd','mt','mtStart','mtEnd','orientation','blastPercentID','genomeLen','mitoLen','assemLen']
    
    nl = '\t'.join(nl) + '\n'
    outFile.write(nl)
    
    assemNum = 1
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        
        extCoords = line[0] + ':' + line[1] + '-' + line[2]
        extFileName = myData['fragAlignDir'] + extCoords + '.fa'
        extBlast = extFileName + '.out'
        extBlastParse = extBlast + '.parse'
        
        assemID = myData['genomeName'] + '_' + 'numtassem' + str(assemNum)
        
        assemLen = int(line[2]) - int(line[1]) + 1
        hspNum = 1
        parseIn = open(extBlastParse,'r')
        for row in parseIn:
            row = row.rstrip()
            row = row.split()
            hspID = assemID + '_hsp' + str(hspNum)
                        
            genomeLen = int(row[2]) - int(row[1]) + 1

            if ',' in row[4]:
                blockStarts = row[4].split(',')
                blockEnds = row[5].split(',')
                blockStarts = [int(j) for j in blockStarts]
                blockEnds = [int(j) for j in blockEnds]                
                mitoLen = (blockEnds[0] - blockStarts[0] + 1) + (blockEnds[1] - blockStarts[1] + 1)                 
            else:
                mitoLen = int(row[5]) - int(row[4]) + 1
            
            nl = [assemID,extCoords,hspID,row[0],row[1],row[2],row[3],row[4],row[5],row[8],row[6],genomeLen,mitoLen,assemLen]
            nl = [str(j) for j in nl]
            nl = '\t'.join(nl) + '\n'
            outFile.write(nl)
            hspNum += 1
        parseIn.close()
        assemNum += 1
    outFile.close()
    inFile.close()
#####################################################################

#### START MAIN PROGRAM ########
parser = argparse.ArgumentParser(description='find numts in assembly')

parser.add_argument('--ref', type=str,help='genome fasta with .fai',required=True)
parser.add_argument('--refname', type=str,help='name of genome',required=True)
parser.add_argument('--mitoref', type=str,help='mitochondria fasta with .fai',required=True)
parser.add_argument('--mitorefdouble', type=str,help='mitochondria fasta concatenated 2x with .fai',required=True)
parser.add_argument('--outdirbase', type=str,help='output directory base dir',required=True)


args = parser.parse_args()
# setup arguments
myData = {}
myData['mitoRefFA'] = args.mitoref
myData['mitoRefFAfai'] = args.mitoref + '.fai'
myData['mitoRefFADouble'] = args.mitorefdouble
myData['mitoRefFADoublefai'] = args.mitorefdouble + '.fai'

myData['genomeName'] = args.refname
myData['genomeFA'] = args.ref
myData['genomeFAfai'] = args.ref + '.fai'
myData['genomeOutDirBase'] = args.outdirbase

# default values
myData['minEval'] = 0.001
myData['chromsToSkipHits'] = ['chrM']
myData['mergeDelta'] = 2000

# check that inputs exist
if os.path.isfile(myData['mitoRefFA']) is False:
    print('ERROR, %s not found!' % myData['mitoRefFA'])
    sys.exit()
if os.path.isfile(myData['mitoRefFAfai']) is False:
    print('ERROR, %s not found!' % myData['mitoRefFAfai'])
    sys.exit()
if os.path.isfile(myData['genomeFA']) is False:
    print('ERROR, %s not found!' % myData['genomeFA'])
    sys.exit()
if os.path.isfile(myData['genomeFAfai']) is False:
    print('ERROR, %s not found!' % myData['genomeFAfai'])
    sys.exit()

if os.path.isfile(myData['mitoRefFADouble']) is False:
    print('ERROR, %s not found!' % myData['mitoRefFADouble'])
    sys.exit()
if os.path.isfile(myData['mitoRefFADoublefai']) is False:
    print('ERROR, %s not found!' % myData['mitoRefFADoublefai'])
    sys.exit()


inFile = open(myData['mitoRefFAfai'],'r')
line = inFile.readline()
line = line.rstrip()
line = line.split()
myData['mitoRefLen'] = int(line[1])    
myData['mitoFAName'] = line[0]
inFile.close()

# check that output dir exists
if os.path.isdir(myData['genomeOutDirBase']) is False:
    print('ERROR %s not found' % myData['genomeOutDirBase'])
    sys.exit()
if myData['genomeOutDirBase'][-1] != '/':
    myData['genomeOutDirBase'] += '/'


# run the analysis    
run_get_genomes_in_assem(myData)
    



