'''
Created on May 21, 2012

@author: lichtens
'''
import csv
import argparse
import sys 
import re
import tempfile
import subprocess 
import string

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

def printStats(individualList):
    for f in individualList:
        r,o = call('cut -f 5,6,11,12 working/' + f + '/' + f + '.maf.oxog.annotated | wc')
        numMutationsPostFiltering = str(int(o.split()[0])-1)
        
        r,o = call('cut -f 5,6,11,12 ~/fhfs/An_SIGMA_Cervical/Individual/' + f + '/jobs/capture/mut/annotated/' + f + '-Tumor.maf.annotated | wc')
        numMutationsPreFiltering = str(int(o.split()[0])-2)
        
        r,o = call('cut -f 5,6,11,12 ~/fhfs/An_SIGMA_Cervical/Individual/' + f + '/jobs/capture/mut/annotated/' + f + '-Tumor.maf.annotated | egrep "\WC\WA" | wc')
        numPossibleCAArtifactMutationsPreFiltering =  str(int(o.split()[0]))
        
        r,o = call('cut -f 5,6,11,12 ~/fhfs/An_SIGMA_Cervical/Individual/' + f + '/jobs/capture/mut/annotated/' + f + '-Tumor.maf.annotated | egrep "\WG\WT" | wc')
        numPossibleGTArtifactMutationsPreFiltering =  str(int(o.split()[0]))
        
        numPossibleArtifactMutationsPreFiltering = str(int(numPossibleCAArtifactMutationsPreFiltering) + int(numPossibleGTArtifactMutationsPreFiltering))
        
        r,o = call('cut -f 5,6,11,12 working/' + f + '/' + f + '.maf.oxog.annotated | egrep "\WC\WA" | wc')
        numPossibleCAArtifactMutationsPostFiltering = o.split()[0]
        
        r,o = call('cut -f 5,6,11,12 working/' + f + '/' + f + '.maf.oxog.annotated | egrep "\WG\WT" | wc')
        numPossibleGTArtifactMutationsPostFiltering = o.split()[0]
        numPossibleArtifactMutationsPostFiltering = str(int(numPossibleCAArtifactMutationsPostFiltering) + int(numPossibleGTArtifactMutationsPostFiltering))
        
        
        print("\t".join([f,numMutationsPreFiltering, numMutationsPostFiltering, numPossibleArtifactMutationsPreFiltering, numPossibleArtifactMutationsPostFiltering]))


if __name__ == '__main__':
    print('WARNING: This is a test script that is likely to fail if not run in lichtens home directory.  Apologies for any inconvenience')
    FAIL = ['CESC-HSCX1127', 'CESC-HSCX1174', 'CESC-HSCX1185', 'CESC-HSCX1215', 'CESC-HSCX1316', 'CESC-HSCX1469', 'CESC-HSCX1484', 'CESC-HSCX1518', 'CESC-HSCX1543', 'CESC-HSCX1550', 'CESC-HSCX1560', 'CESC-HSCX1572', 'CESC-HSCX1598', 'CESC-HSCX1609', 'CESC-HSCX175', 'CESC-HSCX1805', 'CESC-HSCX1856', 'CESC-HSCX2043', 'CESC-HSCX2051', 'CESC-HSCX210', 'CESC-HSCX224', 'CESC-HSCX237', 'CESC-HSCX275', 'CESC-HSCX292', 'CESC-HSCX335', 'CESC-HSCX399', 'CESC-HSCX457', 'CESC-HSCX564', 'CESC-HSCX665', 'CESC-HSCX694', 'CESC-HSCX77', 'CESC-HSCX852', 'CESC-HSCX882', 'CESC-HSCX903', 'CESC-HSCX941', 'CESC-HSCX967', 'CESC-HSCX984']
    OK = ['CESC-HSCX1005', 'CESC-HSCX1025', 'CESC-HSCX1027', 'CESC-HSCX121', 'CESC-HSCX13', 'CESC-HSCX1763', 'CESC-HSCX1777', 'CESC-HSCX1778', 'CESC-HSCX1797', 'CESC-HSCX1989', 'CESC-HSCX2003', 'CESC-HSCX347', 'CESC-HSCX376', 'CESC-HSCX377', 'CESC-HSCX381', 'CESC-HSCX434', 'CESC-HSCX552', 'CESC-HSCX557', 'CESC-HSCX566', 'CESC-HSCX586', 'CESC-HSCX639', 'CESC-HSCX825', 'CESC-HSCX894', 'CESC-HSCX945', 'CESC-HSCX968']
    
    # numMutationsPostFiltering = ' cut -f 5,6,11,12 working/CESC-HSCX1127/CESC-HSCX1127.maf.oxog.annotated | wc'
    # numMutationsPreFiltering = 'cut -f 5,6,11,12 ~/fhfs/An_SIGMA_Cervical/Individual/CESC-HSCX1127/jobs/capture/mut/annotated/CESC-HSCX1127-Tumor.maf.annotated | wc'
    # numPossibleCAArtifactMutationsPreFiltering =  cut -f 5,6,11,12 ~/fhfs/An_SIGMA_Cervical/Individual/CESC-HSCX1127/jobs/capture/mut/annotated/CESC-HSCX1127-Tumor.maf.annotated | egrep "\WC\WA" | wc
    # numPossibleGTArtifactMutationsPreFiltering = cut -f 5,6,11,12 ~/fhfs/An_SIGMA_Cervical/Individual/CESC-HSCX1127/jobs/capture/mut/annotated/CESC-HSCX1127-Tumor.maf.annotated | egrep "\WG\WT" | wc
    # numPossibleCAArtifactMutationsPostFiltering =  cut -f 5,6,11,12 working/CESC-HSCX1127/CESC-HSCX1127.maf.oxog.annotated | egrep "\WC\WA" | wc
    # numPossibleGTArtifactMutationsPostFiltering = cut -f 5,6,11,12 working/CESC-HSCX1127/CESC-HSCX1127.maf.oxog.annotated | egrep "\WG\WT" | wc
    print('ind\tpre\tpost\tpossPre\tpossPost')
    printStats(FAIL)
    
    print('ind\tpre\tpost\tpossPre\tpossPost')
    printStats(OK)
    pass