library("vcfR")

projectname="210813_miseq_AR"
resultsDir="/share/ScratchGeneral/nenbar/projects/Anthony/results"

annotationDir="/share/ScratchGeneral/nenbar/projects/Anthony/annotation"
reference="mhc282.fasta"

vcf_file<-"/share/ScratchGeneral/nenbar/projects/Anthony/results/210813_miseq_AR.vcf"
dna_file <- paste0(annotationDir,"/",reference)

vcf <- read.vcfR( vcf_file, verbose = FALSE )
dna <- ape::read.dna(dna_file, format = "fasta")
#gff <- read.table(gff_file, sep="\t", quote="")

chrom <- create.chromR(name='mhc', vcf=vcf, seq=dna, ann=gff)

plot(chrom)

chrom <- masker(chrom, min_QUAL = 1, min_DP = 300, max_DP = 700, min_MQ = 59.9,  max_MQ = 60.1)
plot(chrom)

chrom <- proc.chromR(chrom, verbose=TRUE)


chromoqc(chrom, dp.alpha=20)
write.vcf()
