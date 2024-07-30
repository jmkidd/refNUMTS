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


#####################################################################
def run_mito_blast2seq(myData):
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
        
#        print(nl)
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
    print('starting merge!')
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
#                        print(tDelt,qDelt)
#                        print(i,hitsPerChrom[tName][i])
#                        print(j,hitsPerChrom[tName][j])
                        numMerges += 1
                        # do merge                        
                        hitsPerChrom[tName][i][2] = hitsPerChrom[tName][j][2]
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
    # and to get the diversity
    print('running align_numt_fragments...')
    myData['fragAlignDir'] = myData['genomeOutDir'] + 'fragalign'
    if os.path.isdir(myData['fragAlignDir'] ) is False:
        print('make frag align dir! %s' % myData['fragAlignDir'] )
        cmd = 'mkdir %s ' % myData['fragAlignDir'] 
        grunCMD(cmd)
    myData['fragAlignDir'] += '/'
    
    inFile = open(myData['blastParseSortMerged'],'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
#        print(line)
        
        extCoords = line[0] + ':' + line[1] + '-' + line[2]
        extFileName = myData['fragAlignDir'] + extCoords + '.fa'
        extBlast = extFileName + '.out'
        extBlastParse = extBlast + '.parse'

        
        cName = line[0]
        cExtStart = int(line[1])
        
        if cExtStart != 100056541:
            continue
        
        segLen = int(line[2]) - int(line[1]) + 1
        
        extCmd = 'samtools faidx %s %s > %s ' % (myData['genomeFA'],extCoords,extFileName)
        runCMD(extCmd)
        
        
        cmd = 'blastn -task blastn -reward 2 -penalty -3 -gapopen 5 -gapextend 2 -evalue 0.001 -outfmt 7 -dust no -query %s -subject %s -out %s \n' % (myData['mitoRefFA'], extFileName,extBlast)
        runCMD(cmd)
        
        hits = read_in_blast_hits(extBlast)
        
        # go through and update them all to have proper coords based on extraction
        for i in range(len(hits)):
            hits[i][0] = cName
            hits[i][1] = hits[i][1] + cExtStart -1
            hits[i][2] = hits[i][2] + cExtStart -1        

        print(len(hits))
        print(hits)
        if len(hits) == 1:
           nl = hits[0]
           nl = [str(i) for i in nl]
           nl = '\t'.join(nl) + '\n'
           outFile = open(extBlastParse,'w')
           outFile.write(nl)
           outFile.close()
        else:
            print(hits)
            sys.exit()

        
        # TO DO NEXT
        # parse hits, get final coordinates, etc...
        # add check that programs are in the proper path, etc..
        
        
    inFile.close()
    


#####################################################################
def read_in_blast_hits(fileName):
    # read in blast hits
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

#### START MAIN PROGRAM ########
parser = argparse.ArgumentParser(description='find numts in assembly')


parser.add_argument('--ref', type=str,help='genome fasta with .fai',required=True)
parser.add_argument('--refname', type=str,help='name of genome',required=True)
parser.add_argument('--mitoref', type=str,help='mitochondria fasta with .fai',required=True)
parser.add_argument('--outdirbase', type=str,help='output directory base dir',required=True)



args = parser.parse_args()
# setup arguments
myData = {}
myData['mitoRefFA'] = args.mitoref
myData['mitoRefFAfai'] = args.ref + '.fai'
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


inFile = open(myData['mitoRefFAfai'],'r')
line = inFile.readline()
line = line.rstrip()
line = line.split()
myData['mitoRefLen'] = int(line[1])    
inFile.close()

# check that output dir exists
if os.path.isdir(myData['genomeOutDirBase']) is False:
    print('ERROR %s not found' % myData['genomeOutDirBase'])
    sys.exit()
if myData['genomeOutDirBase'][-1] != '/':
    myData['genomeOutDirBase'] += '/'


# run the analysis    
run_get_genomes_in_assem(myData)
    


