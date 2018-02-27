#!/bin/bash
# Filter OxoG artifacts from MAF file of called somatic SNVs

# DEPENDENCIES ###############################################################

# NOTE: Edit the below paths according to your environment.

# Matlab (2013a) Common Runtime installation directory
mcr_root=/broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2013a

# Java 1.7
java=/broad/software/free/Linux/redhat_6_x86_64/pkgs/jdk_1.7.0-71/bin/java

# Python 2.7
# Ensure that `python` is is available on $PATH
PATH=/broad/software/free/Linux/redhat_6_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin/:$PATH


# INPUTS #####################################################################

if $# < 6; then
  echo "usage: ${##*/} <id> <bam> <maf> <ref> <dbsnp> <outdir>"
  exit 1
fi

id=$1            # sample id
bam=$2           # path to the BAM file
maf=$3           # path to the MAF file containing called variants
ref=$4           # path to the reference fasta
dbsnp=$5         # path to the dbSNP reference
outdir=$6        # output directory


# PREAMBLE ###################################################################

mkdir -p $outdir

# Directory in which this script resides
root="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


# ANALYSIS ###################################################################

# Task 1
# Calculate the overall OxoG quality score across all reads

# This step takes a long time: skip if output already available
if [[ ! -f "$outdir/${id}.oxoQ.txt" ]]; then
  $root/script/1_annotatePicardOxoQ/annotatePicardOxoQ.sh \
     -J $java \
     -j $root/jar/CollectOxoGMetrics-1.5.jar \
     -i $id \
     -b $bam \
     -c CCG \
     -r $ref \
     -d $dbsnp \
     -o $outdir
fi


# Task 2
# Create interval list at sites of called SNVs

$root/script/2_createOxoGIntervalList/createOxoGIntervals.sh \
  $maf \
  $outdir/${id}.oxoG.interval_list \


# Task 3
# Append the OxoG quality score to the MAF file (inefficient)

value=$(cat "$outdir/${id}.oxoQ.txt")

$root/script/3_AppendAnnotation2MAF/AppendAnnotation2MAF.sh \
   -i ${id} \
   -m ${maf} \
   -f picard_oxoQ \
   -v $value \
   -o $outdir \


# Task 4
# Count reads in the F1R2 and F2R1 configurations
# supporting ref and alt alleles

$java -Xmx2g -jar $root/jar/GenomeAnalysisTK-1.5.jar \
  --analysis_type OxoGMetrics \
  -R ${ref} \
  -I ${bam} \
  -L $outdir/${id}.oxoG.interval_list \
  -o $outdir/${id}.oxoG.counts.txt \


# Task 5
# Append the relevant read counts to the MAF file

$root/script/5_appendOxoGInfo/appendOxoGInfo.sh  \
    --onlyAddColumnsToCopy $outdir/${id}.oxoG.counts.txt \
    $outdir/${id}.picard_oxoQ.maf.annotated \
    $outdir/${id}.oxoG.maf.annotated \


# Task 6
# Apply OxoG filter based on previously collected statistics

inputFile=$outdir/${id}.oxoG.maf.annotated
outputFilename=${id}.oxoG.maf.filtered
poxoG=0.96
artifactThresholdRate=0.01
logThresholdRate=-1
oxoQ_param1=36
oxoQ_param2=1.5

$root/script/6_oxoGFilter_v3/build/run_startFilterMAFFile.sh $mcr_root \
  ${inputFile} ${outputFilename} $outdir \
  0 0 \
  ${poxoG} \
  ${artifactThresholdRate} \
  ${logThresholdRate} \
  ${oxoQ_param1} \
  ${oxoQ_param2} \

# Final output is:
# ${id}.oxoG.maf.filtered

