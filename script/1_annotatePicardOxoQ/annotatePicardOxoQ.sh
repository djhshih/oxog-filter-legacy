#!/bin/sh

# use Python-2.7 
# use Java-1.7

while getopts "i:b:c:r:d:o:j:J:?" Option
do
    case $Option in
        i    ) ID=$OPTARG;;
        b    ) BAM=$OPTARG;;
        c    ) CTX=$OPTARG;;
        r    ) REF=$OPTARG;;
        d    ) DBSNP=$OPTARG;;
        o    ) OUT=$OPTARG;;
        j    ) JAR=$OPTARG;;
        J    ) JAVA=$OPTARG;;
        ?    ) echo "Invalid option: -$OPTARG" >&2
               exit 0;;
        *    ) echo ""
               echo "Unimplemented option chosen."
               exit 0;;
    esac
done


# Have script fail if error is created in subtask.
set -e

echo "ID:       	${ID}"
echo "Bam:         	${BAM}" 
echo "Reference:	${REF}" 
#echo "dbSNP:    	${DBSNP}" 
echo "Context:     	${CTX}" 
echo "Output area:	${OUT}" 
echo "Jar file:	        ${JAR}" 


oopt="-o $OUT "
if [[ -z $OUT ]]; then
   echo "no output full path"
   oopt=""
   OUT="."
fi   

copt="-c $CTX "
if [[ -z $CTX ]]; then
   echo "no context "
   copt=""
fi   
# parse bam path to construct oxog_metric path 
bamdir=$(dirname $BAM)
bamfile=$(basename $BAM)
stub="${bamfile%.*}"
oxoQ=${stub}.oxog_metrics
    
    
Dir=`dirname $0`


if [ ! -f $bamdir/$oxoQ ]; then

    echo "generate oxog_metrics file from bam"
        
    if [[ -n ${DBSNP} && ${DBSNP} != "null" ]]; then
       DB_SNP_ARG="DB_SNP=${DBSNP}"
    else
       DB_SNP_ARG=""
    fi
    ${JAVA} -Xmx3600M -jar ${JAR} INPUT=${BAM} OUTPUT=${OUT}/${oxoQ} REFERENCE_SEQUENCE=${REF} ${DB_SNP_ARG} MINIMUM_QUALITY_SCORE=20 MINIMUM_MAPPING_QUALITY=30 MINIMUM_INSERT_SIZE=60 MAXIMUM_INSERT_SIZE=600 USE_OQ=true CONTEXT_SIZE=1 STOP_AFTER=100000000 VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=LENIENT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false

    if [[ $? -ne 0 ]] ; then
        exit 1
    fi
else
	echo "link to existing oxog_metrics file with bam"
	ln -s  $bamdir/${oxoQ} ${OUT}/.

fi

# oxog_metric file in output area
OXOQ=${OUT}/${oxoQ}

ls -latr ${OXOQ}



echo ""
echo "annotate Picard OxoG  command line: "
echo "python $Dir/annotatePicardOxoQ.py -i $ID -m $OXOQ  $oopt $copt "

python $Dir/annotatePicardOxoQ.py -i $ID -m $OXOQ  $oopt $copt  

echo ""

echo "original report count"

cut -f1-5,11-12 $OXOQ

echo "oxoG: ${CTX}"
cat ${OUT}/${ID}.oxoQ.txt

echo " "
echo "done"
