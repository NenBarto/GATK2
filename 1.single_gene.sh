
conda init 
conda activate gatk

projectname="210813_miseq_AR"

inDir="/share/ScratchGeneral/nenbar/projects/Anthony/raw_data/210813_miseq_AR/fastq/fastq"
resultsDir="/share/ScratchGeneral/nenbar/projects/Anthony/results"
mkdir -p $resultsDir

annotationDir="/share/ScratchGeneral/nenbar/projects/Anthony/annotation"
index=$annotationDir/"mhc282"
logDir="/share/ScratchGeneral/nenbar/projects/Anthony/logs"
mkdir -p $logDir

GATKimage="/share/ScratchGeneral/nenbar/local/images/gatk_latest.sif"

cutadapt_minlen=30

#trimming
revcompDir="$resultsDir/"$projectname".revcomp"
mkdir -p $revcompDir
filterDir="$resultsDir/"$projectname".filter"
mkdir -p $filterDir
mergedDir="$resultsDir/"$projectname".merged"
mkdir -p $mergedDir
cutadaptDir="$resultsDir/"$projectname".cutadapt"
mkdir -p $cutadaptDir
dedupDir="$resultsDir/"$projectname".dedup"
mkdir -p $dedupDir
bwaDir="$resultsDir/"$projectname".bwa"
mkdir -p $bwaDir
haplotypesDir="$resultsDir/"$projectname".haplotypes"
mkdir -p $haplotypesDir
tempRecabBamDir="$resultsDir/"$projectname".tempRecabBam"
mkdir -p $tempRecabBamDir
recabBamDir="$resultsDir/"$projectname".recabBam"
mkdir -p $recabBamDir
sampleDir="$resultsDir/"$projectname".samples"
mkdir -p $sampleDir
tempGVCFdir="$resultsDir/"$projectname".temp-gvcf"
mkdir -p $tempGVCFdir
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
THREADS=5
MINOVLEN=10 #overlap length; default 10
MAXDIFFPCT=100 #overlap cannot have % mismatches higher than this threshold; default 100
MAXEE=2 #Elbrecht 2016 uses maxee=1
MAXEE2=2
MAXDIFFS=10
minuniquesize=10
REFERENCE="mhc282.fasta"

#load files
files=`ls $inDir/*R1_001.fastq.gz`
files=( $files )


for file in ${files[@]}; do 
  echo $file; 
  #file="/share/ScratchGeneral/nenbar/projects/Anthony/raw_data/210813_miseq_AR/fastq/fastq/7111_S29_L001_R1_001.fastq.gz"

  sampleName=`basename $file | sed s/_L001_R1_001.fastq.gz//`
  echo $sampleName

  file1=$inDir/$sampleName"_L001_R1_001.fastq.gz"
  file2=$inDir/$sampleName"_L001_R2_001.fastq.gz"
  fileFilt1=$filterDir/$sampleName"_L001_R1_001_filt.fastq"
  fileFilt2=$filterDir/$sampleName"_L001_R2_001_filt.fastq"
  mergedFile=$mergedDir/$sampleName"_mergedR1R2.fastq"
  fileTrim1=$cutadaptDir/$sampleName"_mergedR1R2_left.fastq"
  fileTrim2=$cutadaptDir/$sampleName"_mergedR1R2_both.fastq"
  fileDedup=$dedupDir/$sampleName"_mergedR1R2_dedup.fastq"
  bamFile="$bwaDir/$sampleName.bam"
  sortedBamExtenstion="$bwaDir/$sampleName.sorted"
  sortedBamFile="$bwaDir/$sampleName.sorted.bam"
  sortedBamFileName="$sampleName.sorted.bam"
  samFile="$bwaDir/$sampleName.sam"
  haplotypeFile=$sampleName.g.vcf.gz

  #filter
  filter_line1="vsearch --threads "$THREADS" --fastq_filter "$file1" --fastq_truncee "$MAXEE" --fastqout $fileFilt1"
  filter_line2="vsearch --threads "$THREADS" --fastq_filter "$file2" --fastq_truncee "$MAXEE" --fastqout $fileFilt2"

  #merge
  merge_line="vsearch --threads $THREADS --fastq_mergepairs $fileFilt1 --reverse $fileFilt2 --fastqout $mergedFile --sizein --sizeout --fastq_allowmergestagger --fastq_maxdiffs $MAXDIFFS --fastq_minovlen $MINOVLEN --fastq_maxdiffpct $MAXDIFFPCT"  

  #trim
  trim_line1="cutadapt -m $cutadapt_minlen --discard-untrimmed -g $FWD -e $e1 -o $fileTrim1 $mergedFile" 
  trim_line2="cutadapt -m $cutadapt_minlen --discard-untrimmed -a $REVrc -e $e2 -o $fileTrim2 $fileTrim1"

  #dedup
  dedup_command="usearch11.0.667_i86linux32 -fastx_uniques $fileTrim2 -fastqout $fileDedup -minuniquesize 10 -sizeout"

  #set up header for read groups, important or breaks GATK 
  #warning, only executable as a job, not in bash
  HEADER=`printf @RG%sID:%s%sSM:%s%sPL:ILLUMINA '\\t' $sampleName '\\t' $sampleName '\\t'`
  mapping_line="bwa mem -R '"${HEADER}"' $index $fileDedup -o $samFile" 
  sam2bam_line="samtools view -S -b $samFile -o $bwaDir/$sampleName.bam"
  sort_line="samtools sort $bwaDir/$sampleName.bam $sortedBamExtenstion"
  index_line="samtools index $sortedBamFile"
  clean_line="rm -f $bamFile & rm -f $samFile"
  haplotypecaller_line="singularity exec --bind $bwaDir:/bwaDir,$haplotypesDir:/haplotypesDir,$annotationDir:/annotationDir $GATKimage gatk --java-options "-Xmx4G" HaplotypeCaller -R /annotationDir/$REFERENCE -I /bwaDir/$sortedBamFileName -O /haplotypesDir/$haplotypeFile -ERC GVCF"

  qsub -b y -wd $logDir -j y -N filter1$sampleName -R y -pe smp 1 -V $filter_line1
  qsub -b y -hold_jid filter1$sampleName -wd $logDir -j y -N filter2$sampleName -R y -pe smp 1 -V $filter_line2
  qsub -b y -hold_jid filter2$sampleName -wd $logDir -j y -N merge$sampleName -R y -pe smp 1 -V $merge_line
  qsub -b y -hold_jid merge$sampleName -wd $logDir -j y -N trim1$sampleName -R y -pe smp 1 -V $trim_line1
  qsub -b y -hold_jid trim1$sampleName -wd $logDir -j y -N trim2$sampleName -R y -pe smp 1 -V $trim_line2
  qsub -b y -hold_jid trim2$sampleName -wd $logDir -j y -N dedup$sampleName -R y -pe smp 1 -V $dedup_command
  qsub -b y -hold_jid dedup$sampleName -wd $logDir -j y -N mapping$sampleName -R y -pe smp 1 -V $mapping_line
  qsub -b y -hold_jid mapping$sampleName -wd $logDir -j y -N sam2bam$sampleName -R y -pe smp 1 -V $sam2bam_line
  qsub -b y -hold_jid sam2bam$sampleName -wd $logDir -j y -N sort$sampleName -R y -pe smp 1 -V $sort_line
  qsub -b y -hold_jid sort$sampleName -wd $logDir -j y -N index$sampleName -R y -pe smp 1 -V $index_line
  qsub -b y -hold_jid index$sampleName -wd $logDir -j y -N clean$sampleName -R y -pe smp 1 -V $clean_line
  qsub -b y -hold_jid clean$sampleName -wd $logDir -j y -N haplotypes$sampleName -R y -pe smp $THREADS -V $haplotypecaller_line

done
