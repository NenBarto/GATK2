
#based on: http://www.ddocent.com/filtering/

#SNPs were removed using VariantFiltration of GATK if: QD <2.0, FS > 40.0, SOR > 5.0, MQ < 20.0, −3.0 > MQRandkSum > 3.0, −3.0 > ReadPosRankSum > 3.0 and AN < 46 (80% of all Alpine ibex individuals). 
#https://www.biorxiv.org/content/10.1101/2020.10.27.357194v1.full


projectname="210813_miseq_AR"
resultsDir="/share/ScratchGeneral/nenbar/projects/Anthony/results/"
annotationDir="/share/ScratchGeneral/nenbar/projects/Anthony/annotation/"
vcfDir="/share/ScratchGeneral/nenbar/projects/Anthony/results/210813_miseq_AR.vcf/"
outDir=$resultsDir/"finalVCF"

mkdir -p $outDir

reference="mhc282.fasta"
inFile=$vcfDir"initial.vcf.gz"


#three step filter
#1. keep variants that have been successfully genotyped in 50% of individuals
#2. a minimum quality score of 30
#3. a minor allele count of 3 - present in at least homo+ het or three hets


vcftools --gzvcf $inFile --max-missing 0.5 --mac 3 --minQ 30 --recode --recode-INFO-all --out $outDir/raw.g5mac3
#12 out of 46 sites

#vcftools --vcf raw.g5mac3.recode.vcf --freq --out clean
vcftools --vcf $outDir/raw.g5mac3.recode.vcf --minDP 10 --recode --recode-INFO-all --out $outDir/raw.g5mac3dp10

#vcftools --vcf raw.g5mac3dp3.recode.vcf --missing-indv
vcftools --vcf $outDir/raw.g5mac3dp10.recode.vcf  --max-missing 0.95 --maf 0.05 --recode --recode-INFO-all --out $outDir/DP10g95maf05 --min-meanDP 20
