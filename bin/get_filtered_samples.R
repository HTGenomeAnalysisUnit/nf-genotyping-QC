#!/bin/env Rscript

library(data.table)
args = commandArgs(trailingOnly=TRUE)

#load data
het <- fread("QC_samples.het.het", header=T, check.names=T)
mis <- fread("QC_samples.miss.imiss", header=T, check.names=T)
sex <- fread("QC_samples.sex_check.sexcheck", header=T, check.names=T)
ibd <- fread("QC_samples.ibd.genome", header=T, check.names=T)

mishet <- data.table(
    FID=het$FID, 
    IID=het$IID, 
    het.rate=(het$N.NM - het$O.HOM)/het$N.NM, 
    mis.rate=mis$F_MISS)

miss.threshold = args[1]
lower.het.threshold = args[2]
upper.het.threshold = args[3]
ibd.threshold = args[4]

write.table(
    mishet[mishet$mis.rate > miss.threshold, .(FID, IID)], 
    file="fail_mis_qc.txt", 
    row.names=F, col.names=F, quote=F)
write.table(
    mishet[mishet$het.rate < lower.het.threshold | mishet$het.rate > upper.het.threshold, .(FID, IID)], 
    file="fail_het_qc.txt", 
    row.names=F, col.names=F, quote=F)
