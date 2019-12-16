load('RData/SKAT_sig.RData')
gene20kb <- read.table('tair10_gene_local_20kb')

geo <- read.table('geo/gene_expression.tsv', header = T, row.names = 1)
tmp <- merge(results_sig, gene20kb, by.x= 2, by.y = 1)
local_eqtl <- tmp[tmp$SetID == tmp$V2,]

## Figure1.R
y <- merge(as.data.frame(table(results_sig$transcript)),
           as.data.frame(table(local_eqtl$transcript)-1),
           by.x =1, by.y=1, all.x = T)
y[is.na(y)] <- 0
x <- y[,2] - y[,3]
x[x>200] <- 210

png('Figure1.png', wid=800, hei=400)
h <- hist(x, main='', ylab='Frequency', xlab='Number of eQTL Genes Identified by One Transcript')
text(h$mids[21], 2400, '>200')
dev.off()


## Figure2.R
gene_loc <- read.table('tair10_gene_28775.bed')

x1 <- merge(results_sig, gene_loc, by.x = 1, by.y = 4)[, 1:4]
colnames(x1)[3:4] <- c('SetID.chrom', 'SetID.start')

x2 <- merge(x1, gene_loc, by.x = 2, by.y = 4)[, 1:6]
colnames(x2)[5:6] <- c('transcript.chrom', 'transcript.start')

chrom_size <- c(30427671, 19698289, 23459830, 18585056, 26975502)
offset <- c(0, cumsum(chrom_size))

for (i in 1:5) {
  x2$SetID.start[x2$SetID.chrom == i] <- x2$SetID.start[x2$SetID.chrom == i] + offset[i]
  x2$transcript.start[x2$transcript.chrom == i] <- x2$transcript.start[x2$transcript.chrom == i] + offset[i]
}

png('Figure2A.png', wid=800, hei=800)
plot(x2$SetID.start, x2$transcript.start, pch=16, cex=0.01, col='blue', xlab='Genomic Location of eQTL Genes', ylab='Genomic Location of Phenotypic Genes', xaxt='n', yaxt='n', bty='n')
abline(h=offset, lwd=2)
abline(v=offset, lwd=2)
mtext(text=paste0('Chr', 1:5), side=1, at=offset[1:5]+diff(offset)/2)
mtext(text=paste0('Chr', 1:5), side=2, at=offset[1:5]+diff(offset)/2)
dev.off()


window_50kb <- read.table('tair10_50kb_windows_genes', stringsAsFactors = F)
x3 <- merge(window_50kb, results_sig, by.x = 1, by.y = 1)
dfm <- as.data.frame(table(x3$V2))

win <- merge(data.frame(win=1:2385), dfm, by.x=1, by.y = 1, all.x = T)
win[is.na(win)] <- 0

window_50kb$chrom <- sub('G.+', 'G', window_50kb$V1)
y <- sapply(split(window_50kb$V2, window_50kb$chrom), function(x){rev(x)[1]})

png('Figure2B.png', wid=800, hei=400)
plot(win$win, win$Freq, type = 'l', xlab='Genomic Location', ylab='Number of Associated Transcripts', col='gray', xaxt='n', las=2)
abline(h=quantile(win$Freq, 0.99), col='red')
abline(v=c(0,y))
mtext(text=paste0('Chr', 1:5), side=1, at=diff(c(0,y))/2 + c(0, y[1:4]), line = 1)
dev.off()

## Figure4.R
mito <- x2[grep('MG', x2$transcript), ]
x4 <- merge(as.data.frame(table(mito$SetID)), gene_loc, by.x = 1, by.y = 4)[, 1:4]
for (i in 1:5) {
  x4$V2[x4$V1 == i] <- x4$V2[x4$V1 == i] + offset[i]
}

png('Figure4A.png', wid=800, hei=400)
plot(x4$V2, x4$Freq, xlim=range(offset), xlab='', ylab='Number of Transcripts', xaxt='n', main=paste0('Mitochondrial Gene (n = ',length(unique(mito$transcript)),')'))
abline(v=offset, lwd=2)
mtext(text=paste0('Chr', 1:5), side=1, at=offset[1:5]+diff(offset)/2, line = 1)
dev.off()

chlo <- x2[grep('CG', x2$transcript), ]
x5 <- merge(as.data.frame(table(chlo$SetID)), gene_loc, by.x = 1, by.y = 4)[, 1:4]
for (i in 1:5) {
  x5$V2[x5$V1 == i] <- x5$V2[x5$V1 == i] + offset[i]
}

png('Figure4B.png', wid=800, hei=400)
plot(x5$V2, x5$Freq, xlim=range(offset), xlab='', ylab='Number of Transcripts', xaxt='n', main=paste0('Chloroplastic Gene (n = ',length(unique(chlo$transcript)),')'))
abline(v=offset, lwd=2)
mtext(text=paste0('Chr', 1:5), side=1, at=offset[1:5]+diff(offset)/2, line = 1)
dev.off()

