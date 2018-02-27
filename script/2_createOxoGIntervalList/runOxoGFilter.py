#!/bin/env python

'''
Created on May 20, 2012

@author: lichtens
'''

import csv
import argparse
import sys 
import re
import tempfile
import subprocess 
import string

WORKSPACE='Lee_NewFilter_Sandbox'
BSUB_PREFIX='bsub -q hour -P  cgafolk -n 4 -R "rusage[mem=4]span[hosts=1]" -o log/'
BSUB_SINGLE_CORE_PREFIX = 'bsub -q hour -P  cgafolk -R "rusage[mem=2]span[hosts=1]" -o log/'

if not (sys.version_info[0] == 2  and sys.version_info[1] in [7]):
    raise "Must use Python 2.7.x"

def call(command, isPrintCmd=False):
    ''' returns returncode, output 
    '''
    try:
        if isPrintCmd:
            print command
        return 0, subprocess.check_output(command, shell=True)
    except subprocess.CalledProcessError as cpe:
        return cpe.returncode,cpe.output

def getIndividualNameFromSampleName(s):
    return re.sub('\-Normal|\-Tumor','', s)

def createIndividualToAnnotatedMafDict(sampleNames):
    result = dict()
    individualNames = set()
    for s in sampleNames:
        individualName = getIndividualNameFromSampleName(s)
        individualNames.add(individualName)
    
    command = 'fiss annot_get ' + WORKSPACE + ' ind=' + ",".join(individualNames) + " maf_file_capture"  
    r,o = call(command)
    
    resultList = o.split('\n')
    for item in resultList:
        tmpList = item.split('\t')
        
        if ( len(tmpList) > 1) and (tmpList[1] <> None ) and (tmpList[1] <> '') and (tmpList[1].startswith('/xchip')):
            result[tmpList[0]] =tmpList[1]
      
    return result
    
if __name__ == '__main__':
    
    # Do fiss calls to setup input files for the rest of the process
    command ='fiss sample_list ' + WORKSPACE + ' PR_SIGMA_Cervical_Capture'
    r, output = call(command)
    
    sampleList = output.split('\n')
    
    # Now that we have the samples, get the bam file for each one.
    command = 'fiss annot_get ' + WORKSPACE + ' samp=' + ",".join(sampleList[1:-1]) + " clean_bam_file_capture "
    r,output = call(command)
    
    
    sampleRows = output.split('\n')
    
    # Get the individuals and find the annotated maf file.
    individualToAnnotatedMAFDict = createIndividualToAnnotatedMafDict(sampleList[1:])
    print individualToAnnotatedMAFDict.keys()
    # Generate the interval list from the annotated maf.
    for ind in individualToAnnotatedMAFDict.keys():
        inputAnnotatedMafFile = individualToAnnotatedMAFDict[ind]
        call('mkdir working/' + ind + "/")
        intervalOutputFilename = 'working/' + ind + "/" + ind + '.interval_list'
        intervalOutputRejectFilename = 'working/'  + ind + "/" + ind +  '.interval_list.reject.txt'
        
        command = 'python readNumPostCallFilter.py ' + inputAnnotatedMafFile + ' ' + intervalOutputFilename + ' ' + intervalOutputRejectFilename
#        call(command)
        print("Already ran: " + command)
        
    # Now that we have the interval list, so we have a way to run OxoGMetrics
    # Create the oxo metrics
    #java -Xmx1g -jar ../../../gatk/dist/GenomeAnalysisTK.jar -nt 4 \
    #   --analysis_type OxoGMetrics -R /seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta -I $pathbam -L interval_File -o outputFile
    
    # For each sample row, run the oxoGMetric gatk task in bsub
    for sampleRow in sampleRows[1:]:
        sampleList = sampleRow.split('\t')
        sampleName = sampleList[0]
        if (len(sampleList) < 2) or (sampleList[1] is None) or sampleList[1] == '':
            continue 
        
        # Determine the interval file to grab
        ind = getIndividualNameFromSampleName(sampleList[0])
        intervalFileName = 'working/' + ind + "/" + ind + '.interval_list'
        
        # Determine the bam file to grab
        bamFile = sampleList[1]
        
        # Create output filename
        outputFilename = 'working/' + ind + "/" + sampleName + '.oxoG.txt'
        
        # Assemble and execute the GATK OxoG command 
        command = BSUB_PREFIX + sampleName + '.log java -Xmx1g -jar ../../dist/GenomeAnalysisTK.jar -nt 4 -T OxoGMetrics -R /seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta -I ' + bamFile + ' -L ' + intervalFileName + ' -o ' + outputFilename
        print(command)
        print('ALREADY DONE')
#        r,o = call(command)
#        print o

    
    for ind in individualToAnnotatedMAFDict.keys():
        
        sampleName = ind + "-Tumor"
        oxoGInputFile = 'working/' + ind + "/" + sampleName + '.oxoG.txt'
        inputFile = individualToAnnotatedMAFDict[ind]    
        
        outputFile = 'working/' + ind + '/' + ind + ".maf.oxog.annotated"
        
        # Assemble and execute command to create the new annotated maf file
        command = BSUB_SINGLE_CORE_PREFIX + sampleName + '.filtering.log python filterByReadConfig.py ' + oxoGInputFile + " " + inputFile + " " + outputFile
        print command
        r,o = call(command)
        print o
    pass