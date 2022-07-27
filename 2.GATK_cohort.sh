
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

sampleDir="$resultsDir/"$projectname".samples"
mkdir -p $sampleDir
tempGVCFdir="$resultsDir/"$projectname".temp-gvcf"
mkdir -p $tempGVCFdir
vcfDir="$resultsDir/"$projectname".vcf"
mkdir -p $vcfDir
resultDir="$resultsDir/"$projectname".result"
mkdir -p $resultDir


REFERENCE="mhc282.fasta"

#load files
files=`ls $inDir/*R1_001.fastq.gz`
files=( $files )


#DBI import
subset="initial"
rm $sampleDir/gvcfs-for-db-import.sample_map
for VCF in $tempGVCFdir/*.gz; do
    filename=$(basename -- "$VCF");
    name="${filename%%.*}";
    file=`basename $VCF`
    VCF="/tempGVCFdir/"$file""
    echo -e "$name\t$VCF" >> $sampleDir/gvcfs-for-db-import.sample_map;
done

#make a database of SNPs   
dbiImport_line="singularity exec --bind $tempGVCFdir:/tempGVCFdir,$sampleDir:/sampleDir $GATKimage gatk --java-options "-Xmx58g" \
    GenomicsDBImport \
    --genomicsdb-workspace-path $subset \
    --tmp-dir tmp \
    --batch-size 20 \
    --sample-name-map /sampleDir/gvcfs-for-db-import.sample_map \
    --L mhc282:1-282"
$dbiImport_line

# call cohort SNPs
genotype_line="singularity exec --bind $vcfDir:/vcfDir,$annotationDir:/annotationDir,$tempGVCFdir:/tempGVCFdir,$sampleDir:/sampleDir $GATKimage gatk --java-options "-Xmx58g" \
    GenotypeGVCFs \
    -R /annotationDir/$REFERENCE \
    -V gendb://initial \
    -O /vcfDir/$subset.vcf.gz"
$genotype_line

#perform vcftools analysis
vcftools --gzvcf $vcfDir/$subset.vcf.gz --indv-freq-burden --out $resultDir/$subset
vcftools --gzvcf $vcfDir/$subset.vcf.gz --freq --out $resultDir/$subset

#vcftools
vcftools \
    --gzvcf $vcfDir/$subset.vcf.gz \
    --minDP 30 \
    --minQ 30 \
    --minGQ 30 \
    --remove-indels --out $resultDir/$subset.filtered


