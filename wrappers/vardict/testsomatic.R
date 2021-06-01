#!/usr/bin/env Rscript

#args <- commandArgs(trailingOnly = TRUE)

dt <- data.table::fread('cat /dev/stdin')
res <- apply(dt[,c(10:13,28:31,8,9,26,27),with=F],1,function(x){
  t1 <- fisher.test(matrix(x[1:4], nrow=2))
  t2 <- fisher.test(matrix(x[5:8], nrow=2))
  t3 <- fisher.test(matrix(c(max(0,x[9]-x[10]), x[10], max(0,x[11]-x[12]), x[12]), nrow=2))
  return(c(round(t1$p.value, 5),round(t1$estimate, 5),round(t2$p.value, 5),round(t2$estimate, 5),round(t3$p.value, 5),round(t3$estimate, 5)))
})

write.table(data.frame(dt[,1:25], res[1,], res[2,], dt[,26:43], res[3,], res[4,], dt[, 44:dim(dt)[2]], res[5,], res[6,]), file = "", quote = F, sep = "\t", eol = "\n", row.names=F, col.names=F)

# d <- read.table( file('stdin'), sep = "\t", header = F, colClasses=c("character", NA, NA, NA, NA, "character", "character", NA, NA, NA, NA, NA, NA, "character", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "character", NA, "character",  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "character", "character", "character", "character"), col.names=c(1:51) )
# 
# if (nrow(d) > 0){
#     pvalues1 <- vector(mode="double", length=dim(d)[1])
#     oddratio1 <- vector(mode="double", length=dim(d)[1])
#     pvalues2 <- vector(mode="double", length=dim(d)[1])
#     oddratio2 <- vector(mode="double", length=dim(d)[1])
#     pvalues <- vector(mode="double", length=dim(d)[1])
#     oddratio <- vector(mode="double", length=dim(d)[1])
# 
#     for( i in 1:dim(d)[1] ) {
# 	h <- fisher.test(matrix(c(d[i,10], d[i,11], d[i,12], d[i,13]), nrow=2))
# 	pvalues1[i] <- round(h$p.value, 5)
# 	oddratio1[i] <- round(h$estimate, 5)
# 	h <- fisher.test(matrix(c(d[i,28], d[i,29], d[i,30], d[i,31]), nrow=2))
# 	pvalues2[i] <- round(h$p.value, 5)
# 	oddratio2[i] <- round(h$estimate, 5)
# 	tref <- if ( d[i,8] - d[i,9] < 0 ) 0 else d[i,8] - d[i,9]
# 	rref <- if ( d[i,26] - d[i,27] < 0 ) 0 else d[i,26] - d[i,27]
# 	h <- fisher.test(matrix(c(tref, d[i,9], rref, d[i,27]), nrow=2))
# 	pvalues[i] <- round(h$p.value, 5)
# 	oddratio[i] <- round(h$estimate, 5)
#     }
# 
#     write.table(data.frame(d[,1:25], pvalues1, oddratio1, d[,26:43], pvalues2, oddratio2, d[, 44:dim(d)[2]], pvalues, oddratio), file = "", quote = F, sep = "\t", eol = "\n", row.names=F, col.names=F)
# }
