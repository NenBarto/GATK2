conda init 
conda activate gatk
module load centos6.10/gi/java/jdk1.8.0_25
#source activate gatk - could be potentially used instead of the containers

#to do
#GATK
#mark duplicates DONE
##Variant Quality Score Recalibration (VQSR) DONE
#integrate GVCF files 
#genotype GVCFs
#R package

projectname="210813_miseq_AR"

inDir="/share/ScratchGeneral/nenbar/projects/Anthony/raw_data/210813_miseq_AR/fastq/fastq"
resultsDir="/share/ScratchGeneral/nenbar/projects/Anthony/results"
mkdir -p $resultsDir

annotationDir="/share/ScratchGeneral/nenbar/projects/Anthony/annotation"
index=$annotationDir/"assembled_contigs"
logDir="/share/ScratchGeneral/nenbar/projects/Anthony/logs"
mkdir -p $logDir

#trimming
trimDir="$resultsDir/"$projectname".trimgalore"
mkdir -p $trimDir
bwaDir="$resultsDir/"$projectname".bwa"
mkdir -p $bwaDir
mrkdupDir="$resultsDir/"$projectname".md-alignments"
mkdir -p $mrkdupDir
sampleDir="$resultsDir/"$projectname".samples"
mkdir -p $sampleDir
cutadapt_minlen=30
tempGVCFdir="$resultsDir/"$projectname".temp-gvcf"
mkdir -p $tempGVCFdir
tempDir="$resultsDir/"$projectname".tmp"
mkdir -p $tempDir
vcfDir="$resultsDir/"$projectname".vcf"
mkdir -p $vcfDir
resultDir="$resultsDir/"$projectname".result"
mkdir -p $resultDir

FWD="CTGGTGATTCCCTGTGACG"
REV="AATGTCTCCACCGCGTTC"
REVrc="GAACGCGGTGGAGACATT"
e1=0.2 #accept 20% error rate in primer sequence
e2=0.2
cutadapt_minlen=30
REFERENCE="mhc282.fasta"
THREADS=5
GATKimage="/share/ScratchGeneral/nenbar/local/images/gatk_latest.sif"

#load files
files=`ls $inDir/711*R1_001.fastq.gz`
files=( $files )



#DBI import
subset="600"

rm $sampleDir/gvcfs-for-db-import.sample_map
for VCF in $tempGVCFdir/$subset*.gz; do
    filename=$(basename -- "$VCF");
    name="${filename%%.*}";
    file=`basename $VCF`
    VCF="/tempGVCFdir/"$file""
    echo -e "$name\t$VCF" >> $sampleDir/gvcfs-for-db-import.sample_map;
done
  
#FIX TEMP!!!  
dbiImport_line="singularity exec --bind $tempGVCFdir:/tempGVCFdir,$sampleDir:/sampleDir $GATKimage gatk --java-options "-Xmx58g" \
    GenomicsDBImport \
    --genomicsdb-workspace-path subset_$subset \
    --tmp-dir tmp \
    --batch-size 20 \
    --sample-name-map /sampleDir/gvcfs-for-db-import.sample_map \
    --L mhc282:1-282"


# now go ahead and run and set batch size equal to cores on node
# also set max RAM a little low, because there is additional overhead
# involved
genotype_line="singularity exec --bind $vcfDir:/vcfDir,$annotationDir:/annotationDir,$tempGVCFdir:/tempGVCFdir,$sampleDir:/sampleDir $GATKimage gatk --java-options "-Xmx58g" \
    GenotypeGVCFs \
    -R annotationDir/$REFERENCE \
    -V gendb://subset_$subset \
    -O /vcfDir/subset_$subset.vcf.gz"


#vcftools
vcftools \
    --gzvcf $vcfDir/subsample_711.vcf.gz \
    --minDP 30 \
    --minQ 30 \
    --minGQ 30 \
    --remove-indels --out $resultDir/subsample_711

vcftools --gzvcf $vcfDir/output.vcf.gz --freq --out $resultDir/output


vcftools --gzvcf $vcfDir/subset_$subset.vcf.gz --indv-freq-burden --out $resultDir/subset_$subset
vcftools --gzvcf $vcfDir/subsample_711.vcf.gz --freq --out $resultDir/subsample_711

#bqsr

#!/bin/bash

# name of output folder
OUTPUT=bqsr-alignments

## DO NOT EDIT BELOW THIS LINE - this comes as input from GNU parallel on STDIN ##
module load jdk/1.8.0_161
source activate gatk

REFERENCE=$1
SITES=$2
INPUT=$3
FILENAME=$(basename -- "$INPUT")
FILENAME_PART="${FILENAME%.*}"
OUT1=$OUTPUT/$FILENAME.recal_data.table
OUT2=$OUTPUT/$FILENAME.bqsr.bam

# ensure that the output directory exists
mkdir -p $OUTPUT
# NOTE - we're assuming single-threaded operation for each BAM file,
# so we set RAM to 4GB each (16 cores, max)
gatk --java-options "-Xmx4G" BaseRecalibrator \
    -I $INPUT \
    -R $REFERENCE \
    --known-sites $SITES \
    -O $OUT1 \
&& gatk --java-options "-Xmx4G" ApplyBQSR \
    -R $REFERENCE \
    -I $INPUT \
    --bqsr-recal-file $OUT1 \
    -O $OUT2
  qsub -b y -hold_jid rev$sampleName -wd $logDir -j y -N gz$sampleName -R y -pe smp 1 -V $gz_command

done
#map files


#IGV