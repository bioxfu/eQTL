library(SKAT)
library(doParallel)
registerDoParallel(core=40)

# Generate SSD file
File.Bed <- 'plink/SNP_candidates.bed'
File.Bim <- 'plink/SNP_candidates.bim'
File.Fam <- 'plink/SNP_candidates.fam'
File.SetID <- 'snp/SNP_candidates.SetID'
File.SSD <- 'plink/SNP_candidates.SSD'
File.Info <- 'plink/SNP_candidates.SSD.info'
Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info)


fam_files <- dir('plink/fam/', full.names = T)
genes <- sub('.fam', '', dir('plink/fam'))

writeLines(c(""), "SKAT_log.txt")
st <- system.time({
  results <- foreach(i=1:length(fam_files), .verbose = F) %dopar% {
    cat(paste("Starting iteration", i, "\n"), file = 'SKAT_log.txt', append = TRUE)
    
    y <- Read_Plink_FAM(fam_files[i], Is.binary = FALSE)$Phenotype
    obj <- SKAT_Null_Model(y ~ 1, out_type = 'C')
    
    # open SSD and Info file and run SKAT
    SSD.INFO <- Open_SSD(File.SSD, File.Info)
    out <- SKAT.SSD.All(SSD.INFO, obj)
    Close_SSD()

    out$results$transcript <- genes[i] 
    out$results
  }
})


save(list = c('st', 'results'), file = 'RData/SKAT.RData')
