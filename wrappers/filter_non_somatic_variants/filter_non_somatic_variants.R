suppressMessages(library(data.table))
suppressMessages(library(openxlsx))
suppressMessages(library(stringi))
suppressMessages(library(fitdistrplus))

Sys.setenv("R_ZIPCMD" = "zip")

run_all <- function(args){
  #read all params as variables
  output_file <- args[1]
  mut_load_output_file <- args[2]
  var_files <- args[3:length(args)]
  
  
  DF <- read_and_filter_var_tabs(var_files)
  
  write.table(DF,file="VARIANTS.tsv",quote = FALSE,sep="\t",col.names = TRUE,row.names = FALSE)
  openxlsx::write.xlsx(c(res_list,res_list_add),file = paste0(per_sample_results_dir,"/",my_sample,".variants.xlsx"))
  
  
}


read_and_filter_var_tabs <- function(var_files){
  #load table
  var_tab <- lapply(var_files,fread)
  names(var_tab) <- gsub(".*/per_sample_final_var_tabs/tsv_formated//?(.*).variants.tsv","\\1",var_files)
  var_tab <- rbindlist(var_tab,use.names = T,idcol = "sample",fill = T)
  
  #add more info from annotation
  var_tab[,alarm := NULL]
  var_tab[COSMIC == ".",COSMIC_hit_count := 0]
  var_tab[COSMIC != ".",COSMIC_hit_count := stri_count_fixed(COSMIC,",") + 1]
  var_tab[,occurence_in_sample := length(unique(sample)), by = var_name]
  var_tab[,c("chr","pos","var") := tstrsplit(var_name,split = "_")]
  var_tab[,c("ref","alt") := tstrsplit(var,split = "/")]
  var_tab[,pos := as.integer(pos)]
  
  #compute distribution params of beta fit for known HET vars 
  het_var_tab <- var_tab[(`1000g_EUR_AF` > 0.2 | gnomAD_NFE > 0.2) & tumor_variant_freq < 0.8]
  
  # puvodne: pridal jsem (het_var_tab$tumor_variant_freq misto tumor_variant_freq
  # het_var_mean_VAF_tab <- het_var_tab[,.(fitdistrplus::fitdist(tumor_variant_freq,distr = "beta")$estimate,type = c("shape1","shape2")), by = sample]
  
  het_var_mean_VAF_tab <- het_var_tab[,.(fitdistrplus::fitdist(het_var_tab$tumor_variant_freq,distr = "beta")$estimate,type = c("shape1","shape2")),by = sample]
  
  het_var_mean_VAF_tab <- dcast.data.table(het_var_mean_VAF_tab,formula = sample ~ type,value.var = "V1")
  
  #add the params to table
  var_tab <- merge(var_tab,het_var_mean_VAF_tab,by = "sample")
  var_tab[,het_var_prob := 1 - abs(pbeta(tumor_variant_freq,shape1,shape2) - 0.5) * 2]
  var_tab[,het_var_prob := max(het_var_prob),by = var_name]
  
  #add some more info
  var_tab[,c("MuTect", "VarScan2", "VarDict", "LoFreq", "Strelka") := lapply(tstrsplit(Called_by,split = ","),function(x) as.logical(as.numeric(x)))]
  var_tab[,called_sum := MuTect + VarScan2 + VarDict + LoFreq + Strelka]
  var_tab[,indel_length := 0]
  var_tab[variant_type == "deletion",indel_length := nchar(ref)]
  var_tab[variant_type == "insertion",indel_length := nchar(alt)]
  var_tab[,tumor_var_reads := tumor_variant_freq * tumor_depth]
  var_tab <- var_tab[`Annotation source` != "."]
  setnames(var_tab,"Annotation source","annot_source")
  
  #filtering step
  #basic germline
  var_tab[,variant_filter := "somatic_OK"]
  var_tab[(`1000g_EUR_AF` > 0.01 | gnomAD_NFE > 0.01) & COSMIC_hit_count < 1,variant_filter := "germline"]
  var_tab[variant_filter != "germline",occurence_in_gene := .N,by = Feature]
  var_tab[tumor_variant_freq > 0.9,variant_filter := "germline"]
  var_tab[occurence_in_sample > 2,variant_filter := "germline"]
  var_tab[(occurence_in_sample > 1 & COSMIC == "." & `md-anderson` == "."),variant_filter := "germline"]
  
  #basic FP variant 
  var_tab[((grepl("MUC|PCDHG|HLA|TRB|IG|LOC",Gene_symbol) & occurence_in_gene > 10) | occurence_in_gene > 100) & variant_filter != "germline",variant_filter := "probable_FP"]
  var_tab[tumor_variant_freq <= 0.02 & variant_filter != "germline",variant_filter := "probable_FP"]
  var_tab[tumor_variant_freq <= 0.05 & (COSMIC_hit_count < 2 | tumor_depth < 100) & variant_filter != "germline",variant_filter := "probable_FP"]
  var_tab[tumor_var_reads <= 5 & variant_filter != "germline",variant_filter := "probable_FP"]
  var_tab[indel_length >= 6 & (indel_length %% 3 != 0 | indel_length >= 20) & variant_filter != "germline",variant_filter := "probable_FP"]
  
  #strict germline
  var_tab <- var_tab[variant_filter == "somatic_OK" & het_var_prob > 0.5,variant_filter := "possible_germline"]
  
  #strict FP variant 
  var_tab[variant_filter == "somatic_OK" & called_sum < 4,variant_filter := "possible_FP"]
  
  
  var_tab[,c("chr","pos","var","ref","alt","mean","sd","occurence_in_gene") := NULL]
  
  setnames(var_tab,"variant_filter","variant_status")
  setcolorder(var_tab,"variant_status")
  
  
  return(var_tab)
  
}



# develop and test
setwd("/home/rj/4TB/GIT/somatic_only_test/")
# args <- c("all_variants.annotated.processed.tsv","mutation_loads.xlsx",
#           "/home/rj/4TB/GIT/somatic_only_test/per_sample_final_var_tabs/tsv_formated/KL2012.variants.tsv",
#           "/home/rj/4TB/GIT/somatic_only_test/per_sample_final_var_tabs/tsv_formated/LN0358.variants.tsv",
#           "/home/rj/4TB/GIT/somatic_only_test/per_sample_final_var_tabs/tsv_formated/MV2008.variants.tsv"
#           )

args <- c("all_variants.annotated.processed.tsv","mutation_loads.xlsx",
          "/home/rj/4TB/TSO_panel_TMB_PATOLOGIE/SNAKEMAKE/per_sample_final_var_tabs/tsv_formated/MS7010.variants.tsv"
)

#run as Rscript
# 
args <- commandArgs(trailingOnly = T)
run_all(args)



