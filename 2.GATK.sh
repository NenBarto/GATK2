conda init 
conda activate enviroDNA
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
cutadapt_minlen=30

tempGVCFdir="$resultsDir/"$projectname".temp-gvcf"
mkdir -p $tempGVCFdir



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
files=`ls $inDir/*R1_001.fastq.gz`
files=( $files )

for file in ${files[@]}; do 
  echo $file; 
  #file="/share/ScratchGeneral/nenbar/projects/Anthony/raw_data/210813_miseq_AR/fastq/fastq/7111_S29_L001_R1_001.fastq.gz"

  sampleName=`basename $file | sed s/_L001_R1_001.fastq.gz//`
  echo $sampleName
  bamFile="$bwaDir/$sampleName.bam"
  samFile="$bwaDir/$sampleName.sam"

  #perform mapping from 1.map.sh 
  ####### notice the change of the name and fix it for the sorted bam
  mappedBamFile="$sampleName.sorted.bam"
  MDfile1=$sampleName.md.bam
  MDfile2=$sampleName.md.fx.bam
  VQSRfile1=$sampleName.g.vcf.gz

  mem="4G"
#  nCores=5

  md_line1="singularity exec --bind $bwaDir:/bwaDir,$mrkdupDir:/mrkdupDir $GATKimage gatk --java-options "-Xmx$mem" MarkDuplicatesSpark --spark-runner LOCAL --input /bwaDir/$mappedBamFile --output /mrkdupDir/$MDfile1 --conf 'spark.executor.cores=$THREADS'"
  md_line2="singularity exec --bind $bwaDir:/bwaDir,$mrkdupDir:/mrkdupDir,$annotationDir:/annotationDir $GATKimage gatk --java-options "-Xmx$mem" SetNmMdAndUqTags --INPUT /mrkdupDir/$MDfile1 --OUTPUT /mrkdupDir/$MDfile2 --REFERENCE_SEQUENCE /annotationDir/$REFERENCE"
  md_line3="singularity exec --bind $bwaDir:/bwaDir,$mrkdupDir:/mrkdupDir $GATKimage gatk --java-options "-Xmx$mem" BuildBamIndex --INPUT /mrkdupDir/$MDfile2"
  VQSR_line1="singularity exec --bind $mrkdupDir:/mrkdupDir,$tempGVCFdir:/tempGVCFdir,$annotationDir:/annotationDir $GATKimage gatk --java-options "-Xmx4G" HaplotypeCaller -R /annotationDir/$REFERENCE -I /mrkdupDir/$MDfile2 -O /tempGVCFdir/$VQSRfile1 -ERC GVCF"
  rm_line="rm $mrkdupDir/$MDfile1 && rm $mrkdupDir/$MDfile1.bai && rm $mrkdupDir/$MDfile1.sbi &&"

  qsub -b y -hold_jid clean$sampleName -wd $logDir -j y -N md1$sampleName -R y -pe smp $THREADS -V $md_line1
  qsub -b y -hold_jid md1$sampleName -wd $logDir -j y -N md2$sampleName -R y -pe smp $THREADS -V $md_line2
  qsub -b y -hold_jid md2$sampleName -wd $logDir -j y -N md3$sampleName -R y -pe smp $THREADS -V $md_line3
  qsub -b y -hold_jid md3$sampleName -wd $logDir -j y -N VQSR$sampleName -R y -pe smp $THREADS -V $VQSR_line1
  qsub -b y -hold_jid VQSR$sampleName -wd $logDir -j y -N clean2$sampleName -R y -pe smp 1 -V $rm_line

done;