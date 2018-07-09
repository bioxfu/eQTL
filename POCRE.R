library(POCRE)
library(data.table)
rna <- read.table('geo/gene_expression.tsv', header = T, row.names = 1)
rna <- t(rna)
rna_pc <- prcomp(rna, center = T, scale. = T)

argv <- commandArgs(T)
input <- argv[1]
# snp <- as.matrix(fread('vcf/filtered_665_accession_snp_maf0.01_gene_Chr1.012'))
# snp_pos <- read.table('vcf/filtered_665_accession_snp_maf0.01_gene_Chr1.012.pos')
snp <- as.matrix(fread(input))
snp_pos <- read.table(paste0(input, '.pos'))

dim(rna)
dim(snp)
dim(snp_pos)

percent <- 0.01
#percent <- 0.05

result <- pocrescreen(rna_pc$x, snp, maxvar=round(ncol(snp)*percent),  maxcmp=5)
snp_pos2 <- snp_pos[result$retSIdx, ]

write.table(snp_pos2, paste0(input, '.pos.select.', percent), quote=F, sep='\t', row.names = F, col.names = F)
