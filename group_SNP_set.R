load('RData/ROCRE.RData')

snp <- unique(do.call(rbind, results))

bed <- data.frame(chr=snp$V1, start=snp$V2-1, end=snp$V2)

bed <- bed[order(bed$chr, bed$start), ]

write.table(bed, 'snp/SNP_candidates.bed', quote = F, sep = '\t', row.names = F, col.names = F)
