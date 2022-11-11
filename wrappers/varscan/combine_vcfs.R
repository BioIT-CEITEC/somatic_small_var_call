library(data.table)

run_all <- function(args){
  vcf1_file <- args[1]
  vcf2_file <- args[2]
  out_file <- args[3]
  
  rvcf1 <- readLines(vcf1_file)
  rvcf2 <- readLines(vcf2_file)
  
  header <- rvcf1[1:grep("^#CHROM",rvcf1)]
  rvcf1 <- rvcf1[(grep("^#CHROM",rvcf1) + 1):length(rvcf1)]
  rvcf2<- rvcf2[(grep("^#CHROM",rvcf2) + 1):length(rvcf2)]
  
  tab1 <- fread(vcf1_file,skip = "#CHROM")
  setnames(tab1,"#CHROM","CHROM")
  tab1[,CHROM := as.character(CHROM)]
  
  tab2 <- fread(vcf2_file,skip = "#CHROM")
  setnames(tab2,"#CHROM","CHROM")
  tab2[,CHROM := as.character(CHROM)]
  
  tab <- rbind(tab1,tab2)
  setkey(tab)

  write(header,file = out_file)
  fwrite(tab,file = out_file, append = TRUE, sep = "\t")
}

# develop and test 2
# args <- character(2)
# args[1] <- "/mnt/ssd/ssd_1/snakemake/CFGenomics/sequencing_results/projects/TP53/2019_Feb/germline_variant_calling/varscan/3655.germline.snps.vcf"
# args[2] <- "/mnt/ssd/ssd_1/snakemake/CFGenomics/sequencing_results/projects/TP53/2019_Feb/germline_variant_calling/varscan/3655.germline.indels.vcf"
# args[3] <- "/mnt/ssd/ssd_1/snakemake/CFGenomics/sequencing_results/projects/TP53/2019_Feb/germline_variant_calling/varscan/3655.germline.vcf"



# run as Rscript
args <- commandArgs(trailingOnly = T)
run_all(args)
