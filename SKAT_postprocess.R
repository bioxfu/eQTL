library(doParallel)
registerDoParallel(core=20)
load('RData/SKAT.RData')

pvalue <- unlist(sapply(results, function(x){x$P.value}))
pvalue_adj <- p.adjust(pvalue, method='BH')

pvalue_cutoff <- max(pvalue[pvalue_adj < 5e-5])

writeLines(c(""), "SKAT_postprocess_log.txt")
st <- system.time({
  results_sig <- foreach(i=1:length(results)) %dopar% {
    cat(paste("Starting iteration", i, "\n"), file = 'SKAT_postprocess_log.txt', append = TRUE)
    x <- results[[i]]
    x <- x[x$P.value < pvalue_cutoff, c('SetID', 'transcript')]
    t(x)
  }
})

results_sig <- as.data.frame(t(do.call(cbind, results_sig)), stringsAsFactors=F)

save(list = c('results_sig'), file = 'RData/SKAT_sig.RData')
