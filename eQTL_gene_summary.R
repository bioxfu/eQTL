library(igraph)
library(GeneOverlap)
library(magrittr)

geo <- read.table('geo/gene_expression.tsv', header = T, row.names = 1)
all_gene <- rownames(geo)

local <- read.table('tair10_gene_local_20kb', stringsAsFactors = F)
local_pair <- paste(local$V1, local$V2)
local2 <- local[local$V1 != local$V2, ]
distant <- setdiff(all_gene, unique(c(local2$V1, local2$V2)))
local3 <- data.frame(from=c(local2$V1, distant), to=c(local2$V2, distant))
g <- graph_from_data_frame(local3)

p_cutoff <- 0.05
marker_cutoff <- 1
dfm <- read.table(paste0('gene/pheno_vs_eQTL_p', p_cutoff, '_marker', marker_cutoff, '.tsv'), header = T, stringsAsFactors = F)
dfm <- dfm[dfm$gene.e %in% all_gene, ]
dfm_local <- dfm[paste(dfm$gene.p, dfm$gene.e) %in% local_pair,]

gene_p <- unique(dfm$gene.p)
gene_p_local <- unique(dfm_local$gene.p)

## summary table
summary_tab <- data.frame(total = c(length(grep('M', all_gene)), length(grep('C', all_gene)), length(grep('[MC]', all_gene, invert = T)), length(all_gene)),
                          with.eQTL = c(length(grep('M', gene_p)), length(grep('C', gene_p)), length(grep('[MC]', gene_p, invert = T)), length(gene_p)),
                          with.eQTL.local = c(length(grep('M', gene_p_local)), length(grep('C', gene_p_local)), length(grep('[MC]', gene_p_local, invert = T)), length(gene_p_local)))

summary_tab$with.eQTL <- paste0(summary_tab$with.eQTL, ' (', round(summary_tab$with.eQTL / summary_tab$total * 100, 1), '%)')
summary_tab$with.eQTL.local <- paste0(summary_tab$with.eQTL.local, ' (', round(summary_tab$with.eQTL.local / summary_tab$total * 100, 1), '%)')
rownames(summary_tab) <- c('Mitochondrion', 'Cholroplast', 'Nuclear genome', 'Total')
colnames(summary_tab) <- c('Total number of genes', 'Number of genes with eQTL', 'Number of genes with local eQTL')
write.table(summary_tab, paste0('results/summary_table_p', p_cutoff, '_marker', marker_cutoff, '.tsv'), sep = '\t', quote = F, col.names = NA)

### Histgram
lst <- split(dfm, dfm$gene.p)
eqtl_gene <- function(x) {
  length(decompose(induced_subgraph(g, x$gene.e)))
}
n <- length(lst)
gene.e_per_gene.p <- c()
pb <- txtProgressBar(min = 0, max = n, style = 3)
for (i in 1:n) {
  gene.e_per_gene.p[i] <- eqtl_gene(lst[[i]])
  setTxtProgressBar(pb, i)
}

pdf('results/Distribution_of_numbers_of_eQTL_identified_for_each_transcript.pdf')
hist(gene.e_per_gene.p, main = '', xlab='Number of eQTL Genes Identifid by One Transcript')
dev.off()

## AgriGO
gene_for_GO <- sort(table(dfm$gene.e), decreasing = T)[1:1000]
write.table(names(gene_for_GO), 'results/eQTL_gene_for_agriGO', quote = F, row.names = F, col.names = F)

## Enrichment analysis
enrichment <- data.frame(matrix(nrow=2, ncol=3))
colnames(enrichment) <- c('paper', 'my', 'overlap')
rownames(enrichment) <- c('flg22_up', 'flg22_down')

uplist <- read.table('doc/flg22_uplist.txt', stringsAsFactors = F)$V1
uplist.eqtl <- read.table('doc/flg22_uplist_eQTL.txt', stringsAsFactors = F)$V1
dfm_uplist <- dfm[dfm$gene.p %in% uplist,]
lst_all <- split(dfm, dfm$gene.e)
lst_uplist <- split(dfm_uplist, dfm_uplist$gene.e)
lst_uplist <- lst_uplist[lapply(lst_uplist, nrow) > 3]

genes <- names(lst_uplist)
pvalues <- c()
pb <- txtProgressBar(min = 0, max = length(genes), style = 3)
for (i in 1:length(genes)) {
  gene <- genes[i]
  eqtl_in_list <- nrow(lst_uplist[[gene]])
  eqtl_not_in_list <- nrow(lst_all[[gene]]) - eqtl_in_list
  not_eqtl_list <- length(uplist) - eqtl_in_list
  not_eqtl_not_in_list <- length(all_gene) - eqtl_not_in_list
  pvalues[i] <- chisq.test(matrix(c(eqtl_in_list, eqtl_not_in_list, not_eqtl_list, not_eqtl_not_in_list), nrow=2))$p.value
  setTxtProgressBar(pb, i)
}

p <- p.adjust(pvalues, method = 'bonferroni')
genes_sig <- genes[p<0.05]

enrichment[1, ] <- c(length(uplist.eqtl), length(genes_sig), length(intersect(genes_sig, uplist.eqtl)))

##
downlist <- read.table('doc/flg22_downlist.txt', stringsAsFactors = F)$V1
downlist.eqtl <- read.table('doc/flg22_downlist_eQTL.txt', stringsAsFactors = F)$V1
dfm_downlist <- dfm[dfm$gene.p %in% downlist,]
lst_all <- split(dfm, dfm$gene.e)
lst_downlist <- split(dfm_downlist, dfm_downlist$gene.e)
lst_downlist <- lst_downlist[lapply(lst_downlist, nrow) > 3]

genes <- names(lst_downlist)
pvalues <- c()
pb <- txtProgressBar(min = 0, max = length(genes), style = 3)
for (i in 1:length(genes)) {
  gene <- genes[i]
  eqtl_in_list <- nrow(lst_downlist[[gene]])
  eqtl_not_in_list <- nrow(lst_all[[gene]]) - eqtl_in_list
  not_eqtl_list <- length(downlist) - eqtl_in_list
  not_eqtl_not_in_list <- length(all_gene) - eqtl_not_in_list
  pvalues[i] <- chisq.test(matrix(c(eqtl_in_list, eqtl_not_in_list, not_eqtl_list, not_eqtl_not_in_list), nrow=2))$p.value
  setTxtProgressBar(pb, i)
}

p <- p.adjust(pvalues, method = 'bonferroni')
genes_sig <- genes[p<0.05]
enrichment[2, ] <- c(length(downlist.eqtl), length(genes_sig), length(intersect(genes_sig, downlist.eqtl)))
write.table(enrichment, 'results/flg22_enrichment_eQTL_gene_summary.tsv', sep = '\t', quote = F, col.names = NA)


