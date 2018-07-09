geo <- read.table('geo/gene_expression.tsv', header = T, row.names = 1)
pheno <- rownames(geo)

n <- nrow(geo)
#n <- 100
p_cutoff <- 0.05
marker_cutoff <- 1

mat <- matrix(nrow = 2, ncol=n*n)
rownames(mat) <- c('gene.p', 'gene.e')
start <- 1
pb <- txtProgressBar(min = 0, max = n, style = 3)
for (i in 1:n) {
  test <- read.table(paste0('skat/test_', i, '.tsv'), header = T, stringsAsFactors = F)
  p <- p.adjust(test$P.value, method = 'BH')
  marker <- test$N.Marker.Test
  index <- p < p_cutoff & marker >= marker_cutoff
  index[is.na(index)] <- FALSE
  l <- sum(index)
  if (l > 0) {
    gene.e <- test$SetID[index]
    gene.p <- rep(pheno[i], length(gene.e))
    end <- start + l - 1
    mat[, start:end] <- rbind(gene.p, gene.e)
    start <- end + 1
  }
  setTxtProgressBar(pb, i)
}
dfm <- as.data.frame(t(mat[, 1:end]))

write.table(dfm, paste0('gene/pheno_vs_eQTL_p', p_cutoff, '_marker', marker_cutoff, '.tsv'), row.names = F, quote = F, sep = '\t')
