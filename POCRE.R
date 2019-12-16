library(POCRE)
library(data.table)
library(doParallel)
registerDoParallel(core=40)

rna <- read.table('geo/gene_expression.tsv', header = T, row.names = 1)
rna <- t(rna)

snp <- as.matrix(fread('vcf/filtered_665_accession_snp_maf0.01_gene_matrix.012'))
snp_pos <- read.table('vcf/filtered_665_accession_snp_maf0.01_gene_matrix.012.pos')

dim(rna)
dim(snp)
dim(snp_pos)

writeLines(c(""), "POCRE_log.txt")
st <- system.time({
  results <- foreach(i=1:ncol(rna)) %dopar% {
    cat(paste("Starting iteration", i, "\n"), file = 'POCRE_log.txt', append = TRUE)
    result <- pocrescreen(rna[, i, drop=F], snp)
    snp_pos[result$retSIdx, ]
  }
})

save(list = c('st', 'results'), file = 'RData/ROCRE.RData')

