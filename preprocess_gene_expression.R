geo <- read.table('geo/GSE80744_ath1001_tx_norm_2016-04-21-UQ_gNorm_normCounts_k4.tsv', header = T, row.names = 1)
colnames(geo) <- sub('X', '', colnames(geo))
acc <- read.table('665_accession')$V1

geo2 <- round(log2(geo[, colnames(geo) %in% c('gene_id', acc)] + 1), 2)
write.table(geo2, 'geo/gene_expression.tsv', quote=F, sep='\t', col.names = NA)

geo_pheno <- data.frame(FID=colnames(geo2),
                        IID=colnames(geo2),
                        t(geo2))
write.table(geo_pheno, 'geo/gene_expression_pheno.tsv', quote=F, sep='\t', row.names = FALSE)
