'''
Created on May 18, 2012

@author: lichtens
'''

import csv
import argparse
import sys 
import re
import tempfile
import subprocess 

if not (sys.version_info[0] == 2  and sys.version_info[1] in [7]):
    raise "Must use Python 2.7.x"

def call(command, isPrintCmd=True):
    ''' returns returncode, output 
    '''
    try:
        if isPrintCmd:
            print command
        return 0, subprocess.check_output(command, shell=True)
    except subprocess.CalledProcessError as cpe:
        return cpe.returncode,cpe.output

def parseOptions():
    description = '''Given a tsv call_stats file with headers on the first line, run each mutation through the .'''
    epilog= '''
        '''
    parser = argparse.ArgumentParser(description=description, epilog=epilog, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('inputFile', metavar='inputFile',  type=str, help ='call stats file from MuTect or an annotated maf from Oncotator ')
    parser.add_argument('outputIntervalFileName', metavar='outputIntervalFileName', type=str, help ='')
    parser.add_argument('outputIntervalRejectFileName', metavar='outputIntervalRejectFileName', type=str, help ='Only used if input is a call_stats file (or other file with judgement column)')
    
    args = parser.parse_args()
    
    return args 

if __name__ == '__main__':
    
    args = parseOptions()
    inputFile = args.inputFile
    outputIntervalFileName = args.outputIntervalFileName
    outputIntervalRejectFileName = args.outputIntervalRejectFileName
    
    outFP = file(outputIntervalFileName,'w')
    outRejectFP = file(outputIntervalRejectFileName,'w')
    inputFP = file(inputFile,'r')
    inputFP.readline()
    
    inputTSVReader = csv.DictReader(inputFP, delimiter='\t')
       
    # Read in the call stats file and for each chr, start row, add it to the interval file (of length 1)
    #contig  position
    for line in inputTSVReader:
        if 'contig' in line.keys():
            chrom = line['contig']
            position = line['position']
        else:
            chrom = line['Chromosome']
            position = line['Start_position']
                        
        if ('judgement' not in line.keys()) or (line['judgement'] == 'KEEP'):
            outFP.write(chrom + ":" + position + "\n")
        else:
            outRejectFP.write(chrom + ":" + position + "\n")
    
    outFP.close()
    outRejectFP.close()
    pass