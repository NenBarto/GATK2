conda init 
conda activate enviroDNA

#mapping to OTU3 as the most common OTU of length 282, named mhc282
#bwa index -p mhc282 -a is mhc282.fasta
#reference should be created with gatk
#gatk-launch CreateSequenceDictionary -R ref.fasta
#singularity exec --bind $annotationDir:/annotationDir $GATKimage gatk CreateSequenceDictionary -R /annotationDir/$REFERENCE


#to do
#GATK
#R package

projectname="210813_miseq_AR"

inDir="/share/ScratchGeneral/nenbar/projects/Anthony/raw_data/210813_miseq_AR/fastq/fastq"
resultsDir="/share/ScratchGeneral/nenbar/projects/Anthony/results"
mkdir -p $resultsDir

annotationDir="/share/ScratchGeneral/nenbar/projects/Anthony/annotation"
index=$annotationDir/"mhc282"
logDir="/share/ScratchGeneral/nenbar/projects/Anthony/logs"
mkdir -p $logDir

cutadapt_minlen=30

#trimming
revcompDir="$resultsDir/"$projectname".revcomp"
mkdir -p $revcompDir

trimDir="$resultsDir/"$projectname".trimgalore"
mkdir -p $trimDir
bwaDir="$resultsDir/"$projectname".bwa"
mkdir -p $bwaDir

FWD="CTGGTGATTCCCTGTGACG"
REV="AATGTCTCCACCGCGTTC"
REVrc="GAACGCGGTGGAGACATT"
e1=0.2 #accept 20% error rate in primer sequence
e2=0.2
cutadapt_minlen=30

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
  #file2rev=$revcompDir/$sampleName"_L001_R2_001.revcomp.fastq"
  fileTrim1=$trimDir/$sampleName"_L001_R1_001_val_1.fq"
  fileTrim2=$trimDir/$sampleName"_L001_R2_001_val_2.fq"
  fileTrim1gz=$trimDir/$sampleName"_L001_R1_001_val_1.fq.gz"
  fileTrim2gz=$trimDir/$sampleName"_L001_R2_001_val_2.fq.gz"
  bamFile="$bwaDir/$sampleName.bam"
  sortedBamFile="$bwaDir/$sampleName.sorted"
  samFile="$bwaDir/$sampleName.sam"


  #rev_command="vsearch --fastx_revcomp $file2 --fastqout $file2rev"
  trim_line="cutadapt -g $FWD -G $REV --discard-untrimmed -o $fileTrim1 -p $fileTrim2 $file1 $file2"
  gz_command1="gzip -f $fileTrim1"
  gz_command2="gzip -f $fileTrim2"

  #set up header for read groups, important or breaks GATK
  HEADER=`printf @RG%sID:%s%sSM:%s%sPL:ILLUMINA '\\t' $sampleName '\\t' $sampleName '\\t'`
  mapping_line="bwa mem -R '"${HEADER}"' $index $fileTrim1gz $fileTrim2gz -o $samFile" 
  sam2bam_line="samtools view -S -b $samFile -o $bwaDir/$sampleName.bam"
  sort_line="samtools sort $bwaDir/$sampleName.bam $sortedBamFile"
  index_line="samtools index $bwaDir/$sampleName.sorted.bam"
  clean_line="rm -f $bamFile & rm -f $samFile"

  qsub -b y -wd $logDir -j y -N trim$sampleName -R y -pe smp 1 -V $trim_line
  qsub -b y -hold_jid trim$sampleName -wd $logDir -j y -N gz1$sampleName -R y -pe smp 1 -V $gz_command1
  qsub -b y -hold_jid gz1$sampleName -wd $logDir -j y -N gz2$sampleName -R y -pe smp 1 -V $gz_command2
  qsub -b y -hold_jid gz2$sampleName -wd $logDir -j y -N map$sampleName -R y -pe smp 1 -V $mapping_line
  qsub -b y -hold_jid map$sampleName -wd $logDir -j y -N sam2bam$sampleName -R y -pe smp 1 -V $sam2bam_line
  qsub -b y -hold_jid sam2bam$sampleName -wd $logDir -j y -N sort$sampleName -R y -pe smp 1 -V $sort_line
  qsub -b y -hold_jid sort$sampleName -wd $logDir -j y -N index$sampleName -R y -pe smp 1 -V $index_line
  qsub -b y -hold_jid index$sampleName -wd $logDir -j y -N clean$sampleName -R y -pe smp 1 -V $clean_line

done
#map files


#IGV