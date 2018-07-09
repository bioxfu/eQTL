library(SKAT)

argv <- commandArgs(T)
n <- argv[1]
#p <- 0.05
p <- 0.01

cmd <- paste0('./src/plink --memory 100 --vcf vcf/filtered_665_accession_snp_maf0.01_gene.pos.select.', p, '.ensembl.vcf.gz --pheno geo/gene_expression_pheno.tsv --mpheno ', n, ' --make-bed --out plink/test_', n)
system(cmd)

# Generate SSD file
File.Bed <- paste0('plink/test_', n, '.bed')
File.Bim <- paste0('plink/test_', n, '.bim')
File.Fam <- paste0('plink/test_', n, '.fam')
File.SetID <- paste0('tair10_gene.SetID.', p)
File.SSD <- paste0('plink/test_', n, '.SSD')
File.Info <- paste0('plink/test_', n, '.SSD.info')

Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info)

# open SSD and Info file and run SKAT
FAM <- Read_Plink_FAM(File.Fam, Is.binary = FALSE)

y <- FAM$Phenotype

SSD.INFO <- Open_SSD(File.SSD, File.Info)

obj <- SKAT_Null_Model(y ~ 1, out_type = 'C')

out <- SKAT.SSD.All(SSD.INFO, obj)

write.table(out$results, paste0('skat/test_', n, '.tsv'), row.names = F, sep = '\t', quote = F)

Close_SSD()

system(paste0('rm plink/test_', n, '.*'))
