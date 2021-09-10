library(data.table)
library(vcfR)

get_depth_and_bias_from_DP4_format <- function(DP4){
  a <- as.matrix(as.data.table(lapply(tstrsplit(DP4,","),as.numeric)))
  return(list(rowSums(a),(a[,1] + a[,3]) / rowSums(a)))
}

run_all <- function(args){
  snp_var_file <- args[1]
  indel_var_file <- args[2]
  output_file <- args[3]
  calling_type <- args[4]
  
  if(calling_type == "paired"){
    # SNV processing
    vcf <- vcfR::read.vcfR(snp_var_file,verbose = F)
    snp_var_tab <- data.table(chrom = vcf@fix[,"CHROM"]
                              ,position = as.integer(vcf@fix[,"POS"])
                              ,reference = vcf@fix[,"REF"]
                              ,alternative = vcf@fix[,"ALT"]
                              ,filter = vcf@fix[,"FILTER"]
                              ,tumor_variant_freq = vcfR::extract.gt(vcf,element = "VAF",as.numeric = T)[,2]
                              ,tumor_depth = vcfR::extract.gt(vcf,element = "DP4")[,2]
                              ,normal_variant_freq = vcfR::extract.gt(vcf,element = "VAF",as.numeric = T)[,1]
                              ,normal_depth = vcfR::extract.gt(vcf,element = "DP4")[,1]
                              ,number_of_callers = vcfR::extract.info(vcf,"NUM_TOOLS")
                              ,callers = vcfR::extract.info(vcf,"MVSDULK"))
    
    if(nrow(snp_var_tab)){
      snp_var_tab[,c("tumor_depth","tumor_fwd_strand_pct") := get_depth_and_bias_from_DP4_format(tumor_depth)]
      snp_var_tab[,c("normal_depth","normal_fwd_strand_pct") := get_depth_and_bias_from_DP4_format(normal_depth)]
    }
    
    # INDEL processing
    vcf <- vcfR::read.vcfR(indel_var_file,verbose = F)
    indel_var_tab <- data.table(chrom = vcf@fix[,"CHROM"]
                                ,position = as.integer(vcf@fix[,"POS"])
                                ,reference = vcf@fix[,"REF"]
                                ,alternative = vcf@fix[,"ALT"]
                                ,filter = vcf@fix[,"FILTER"]
                                ,tumor_variant_freq = vcfR::extract.gt(vcf,element = "VAF",as.numeric = T)[,2]
                                ,tumor_depth = vcfR::extract.gt(vcf,element = "DP4")[,2]
                                ,normal_variant_freq = vcfR::extract.gt(vcf,element = "VAF",as.numeric = T)[,1]
                                ,normal_depth = vcfR::extract.gt(vcf,element = "DP4")[,1]
                                ,number_of_callers = vcfR::extract.info(vcf,"NUM_TOOLS")
                                ,callers = vcfR::extract.info(vcf,"MVDLK"))
    
    if(nrow(indel_var_tab)){
      indel_var_tab[,c("tumor_depth","tumor_fwd_strand_pct") := get_depth_and_bias_from_DP4_format(tumor_depth)]
      indel_var_tab[,c("normal_depth","normal_fwd_strand_pct") := get_depth_and_bias_from_DP4_format(normal_depth)]
    }
    
    
    
  } else {
    # SNV processing
    vcf <- vcfR::read.vcfR(snp_var_file,verbose = F)
    snp_var_tab <- data.table(chrom = vcf@fix[,"CHROM"]
                              ,position = as.integer(vcf@fix[,"POS"])
                              ,reference = vcf@fix[,"REF"]
                              ,alternative = vcf@fix[,"ALT"]
                              ,filter = vcf@fix[,"FILTER"]
                              ,tumor_variant_freq = vcfR::extract.gt(vcf,element = "VAF",as.numeric = T)[,1]
                              ,tumor_depth = vcfR::extract.gt(vcf,element = "DP4")[,1]
                              ,number_of_callers = vcfR::extract.info(vcf,"NUM_TOOLS")
                              ,callers = vcfR::extract.info(vcf,"MVDLK"))
    
    if(nrow(snp_var_tab)){
      snp_var_tab[,c("tumor_depth","tumor_fwd_strand_pct") := get_depth_and_bias_from_DP4_format(tumor_depth)]
    }
    
    # INDEL processing
    vcf <- vcfR::read.vcfR(indel_var_file,verbose = F)
    indel_var_tab <- data.table(chrom = vcf@fix[,"CHROM"]
                                ,position = as.integer(vcf@fix[,"POS"])
                                ,reference = vcf@fix[,"REF"]
                                ,alternative = vcf@fix[,"ALT"]
                                ,filter = vcf@fix[,"FILTER"]
                                ,tumor_variant_freq = vcfR::extract.gt(vcf,element = "VAF",as.numeric = T)[,1]
                                ,tumor_depth = vcfR::extract.gt(vcf,element = "DP4")[,1]
                                ,number_of_callers = vcfR::extract.info(vcf,"NUM_TOOLS")
                                ,callers = vcfR::extract.info(vcf,"MVDLK"))
    
    if(nrow(indel_var_tab)){
      indel_var_tab[,c("tumor_depth","tumor_fwd_strand_pct") := get_depth_and_bias_from_DP4_format(tumor_depth)]
    }
  }
  
  
  if(nrow(indel_var_tab)){
    indel_var_tab[,new_ref := reference]
    indel_var_tab[sapply(seq_along(reference),function(x) grepl(reference[x],alternative[x])),new_ref := "-"]
    indel_var_tab[sapply(seq_along(reference),function(x) grepl(alternative[x],reference[x])),new_ref := stringi::stri_sub(new_ref,from = nchar(alternative) + 1)]
    
    indel_var_tab[,new_alt := alternative]
    indel_var_tab[sapply(seq_along(reference),function(x) grepl(alternative[x],reference[x])),new_alt := "-"]
    indel_var_tab[sapply(seq_along(reference),function(x) grepl(reference[x],alternative[x])),new_alt := stringi::stri_sub(new_alt,from = nchar(reference) + 1)]
    
    indel_var_tab[,new_pos := position]
    indel_var_tab[sapply(seq_along(reference),function(x) grepl(reference[x],alternative[x])),new_pos := new_pos + nchar(reference)]
    indel_var_tab[sapply(seq_along(reference),function(x) grepl(alternative[x],reference[x])),new_pos := new_pos + nchar(alternative)]
    indel_var_tab[,c("position","reference","alternative") := .(new_pos,new_ref,new_alt)]
    indel_var_tab[,c("new_pos","new_ref","new_alt") := NULL]
    
    var_tab <- rbind(snp_var_tab,indel_var_tab)
  }else{var_tab <- snp_var_tab}
  
  
  # var_tab <- rbind(snp_var_tab,indel_var_tab)
  
  
  
  # PRINT CALLING STATISTICS
  stat_tab <- copy(var_tab)
  
  # if(nrow(var_tab) > 0){
    
    # # GET CALLING STATISTICS
    # stat_tab[,is_pass := filter == "PASS"]
    # stat_tab[!(is.na(status) | status == "LikelySomatic" | status == "StrongSomatic"),status := "NotSomatic"]
    # stat_tab[,status := factor(status,levels = c("NotSomatic","LikelySomatic","StrongSomatic"))]
    # stat_tab[,by_caller_sum := .N,by = c("caller")]
    # 
    # stat_tab_out1 <- stat_tab[,list(variants = .N),by = c("caller")]
    # stat_tab_out1[,c("is_pass","status","percent_of_caller") := NA]
    # stat_tab_out2 <- stat_tab[,list(variants = .N,percent_of_caller = round(.N / by_caller_sum[1] * 100,1)),by = c("caller","is_pass","status")]
    # setcolorder(stat_tab_out1,names(stat_tab_out2))
    # setorder(stat_tab_out1)
    # setorder(stat_tab_out2)
    # 
    # write.table(rbind(stat_tab_out1,stat_tab_out2),file = gsub(".tsv$",".stats_tab.tsv",output_file),sep = "\t",row.names = F,col.names = T,quote = F,na = "")
    # 
    # caller_types <- unique(stat_tab$caller)
    # venn_list <- lapply(caller_types,function(x) stat_tab[caller == x,paste(chrom,position,sep = "_")])
    # is_ok <- venn.diagram(venn_list,gsub(".tsv$",".stats_ven_all.tiff",output_file)
    #                       ,imagetype = "tiff"
    #                       ,fill=c("red", "green","blue","yellow")[seq_along(caller_types)]
    #                       ,alpha=c(0.5,0.5,0.5,0.5)[seq_along(caller_types)]
    #                       ,cex = 1.5
    #                       ,cat.fontface=2
    #                       ,category.names=caller_types,main = "Somatic variants all")
    # 
    # stat_tab <- copy(var_tab)
    # 
    # venn_list <- lapply(caller_types,function(x) stat_tab[caller == x,paste(chrom,position,sep = "_")])
    # is_ok <- venn.diagram(venn_list,gsub(".tsv$",".stats_ven_pass.tiff",output_file)
    #                       ,imagetype = "tiff"
    #                       ,fill=c("red", "green","blue","yellow")[seq_along(caller_types)]
    #                       ,alpha=c(0.5,0.5,0.5,0.5)[seq_along(caller_types)]
    #                       ,cex = 1.5
    #                       ,cat.fontface=2
    #                       ,category.names=caller_types,main = "Somatic variants PASS")
    # 
    # var_tab[,status := gsub("(.).*","\\1",status)]
    # var_tab[caller == "vardict" ,caller := paste0(caller,"_",status)]
  # }
  
  #FILTER CALLS
  var_tab <- var_tab[filter == "PASS"]
  
  if(calling_type == "paired"){
    var_tab <- var_tab[tumor_depth > 0 & tumor_variant_freq > 0 & normal_depth > 0]
  } else {
    var_tab <- var_tab[tumor_depth > 0 & tumor_variant_freq > 0]
  }
  var_tab[tumor_variant_freq > 1,tumor_variant_freq := 1]
  
  var_tab[,filter := NULL]
  
  var_tab[,var_name := paste0(chrom,"_",position,"_",reference,"/",alternative)]
  setorder(var_tab,chrom,position,reference)
  if(calling_type == "paired"){
    var_tab <- var_tab[,.(var_name,tumor_variant_freq,tumor_depth,normal_variant_freq,normal_depth,number_of_callers,callers,tumor_fwd_strand_pct,normal_fwd_strand_pct)]
  } else {
    var_tab <- var_tab[,.(var_name,tumor_variant_freq,tumor_depth,number_of_callers,callers,tumor_fwd_strand_pct)]
  }
  
  fwrite(var_tab,file = output_file,sep = "\t")
  
  
}


# develop and test
#
# args <- c("somatic_seq_results/Pca_25/Consensus.sSNV.vcf","somatic_seq_results/Pca_25/Consensus.sINDEL.vcf","somatic_seq_results/Pca_25.variants.tsv")
# setwd("/mnt/ssd/ssd_1/snakemake/stage359_PC.seq_A/somatic_variant_calling")
# args <- c("somatic_seq_results/PCA_116/Consensus.sSNV.vcf","somatic_seq_results/PCA_116/Consensus.sINDEL.vcf","somatic_seq_results/PCA_116.variants.tsv","paired","")
# setwd("/mnt/ssd/ssd_1/snakemake/stage447_ACGT03A.Pca_E48/somatic_seq")

#run as Rscript
#
args <- commandArgs(trailingOnly = T)
run_all(args)



