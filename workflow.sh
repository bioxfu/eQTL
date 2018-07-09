source activate gmatic

# http://www.stat.purdue.edu/~zhangdb/packages/POCRE/
mkdir src
mkdir doc
mkdir data
wget -c http://www.stat.purdue.edu/~zhangdb/packages/POCRE/POCRE_0.1.0.tar.gz -P src
wget -c http://www.stat.purdue.edu/~zhangdb/packages/POCRE/pocreworkshop.pdf -P doc
wget -c http://www.stat.purdue.edu/~zhangdb/packages/POCRE/metabdata.csv -P data
wget -c http://www.stat.purdue.edu/~zhangdb/packages/POCRE/pe.xlsx -P data
wget -c http://www.stat.purdue.edu/~zhangdb/packages/POCRE/GSE58210_NormalizedData_withannotations.txt -P data

# VCF (http://1001genomes.org/data/GMI-MPI/releases/v3.1/)
mkdir vcf
wget http://1001genomes.org/data/GMI-MPI/releases/v3.1/intersection_snp_short_indel_vcf/ -P vcf
grep 'vcf.gz' vcf/index.html|sed -r 's/.+"intersection_//'|sed -r 's/\.vcf.+//'|sort -n > vcf/vcf_accession_list

# GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80744)
mkdir geo
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE80nnn/GSE80744/suppl/GSE80744_ath1001_tx_norm_2016-04-21-UQ_gNorm_normCounts_k4.tsv.gz -P geo
zcat geo/GSE80744_ath1001_tx_norm_2016-04-21-UQ_gNorm_normCounts_k4.tsv.gz |head -1|sed 's/X//g'|sed -r 's/\t/\n/g'|grep -v 'gene' > geo/geo_accessions_list
zcat geo/GSE80744_ath1001_tx_norm_2016-04-21-UQ_gNorm_normCounts_k4.tsv.gz > geo/GSE80744_ath1001_tx_norm_2016-04-21-UQ_gNorm_normCounts_k4.tsv

# find 665 accession with both VCF and GEO data
cat vcf/vcf_accession_list geo/geo_accessions_list|sort -n|uniq -d > 665_accession

# extract 655 accession in GEO data
Rscript preprocess_gene_expression.R

# extract SNPs (MAF>0.01) within genes in 655 accession from VCF file
wget -c http://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz -P vcf
bcftools index -t vcf/1001genomes_snp-short-indel_only_ACGTN.vcf.gz
bcftools index --nrecords vcf/1001genomes_snp-short-indel_only_ACGTN.vcf.gz # get total variant count
bcftools index --stats vcf/1001genomes_snp-short-indel_only_ACGTN.vcf.gz # get variant count per chromsome

vcftools --gzvcf vcf/1001genomes_snp-short-indel_only_ACGTN.vcf.gz --keep 665_accession --remove-indels --maf 0.01 --bed tair10_gene.bed --recode --recode-INFO-all --stdout | bgzip -c > vcf/filtered_665_accession_snp_maf0.01_gene.vcf.gz
bcftools index -t vcf/filtered_665_accession_snp_maf0.01_gene.vcf.gz
bcftools index --nrecords vcf/filtered_665_accession_snp_maf0.01_gene.vcf.gz # get total variant count
bcftools index --stats vcf/filtered_665_accession_snp_maf0.01_gene.vcf.gz # get variant count per chromsome

# split into five chrom
vcftools --gzvcf vcf/filtered_665_accession_snp_maf0.01_gene.vcf.gz --chr 1 --recode --recode-INFO-all --stdout | bgzip -c > vcf/filtered_665_accession_snp_maf0.01_gene_Chr1.vcf.gz &
vcftools --gzvcf vcf/filtered_665_accession_snp_maf0.01_gene.vcf.gz --chr 2 --recode --recode-INFO-all --stdout | bgzip -c > vcf/filtered_665_accession_snp_maf0.01_gene_Chr2.vcf.gz &
vcftools --gzvcf vcf/filtered_665_accession_snp_maf0.01_gene.vcf.gz --chr 3 --recode --recode-INFO-all --stdout | bgzip -c > vcf/filtered_665_accession_snp_maf0.01_gene_Chr3.vcf.gz &
vcftools --gzvcf vcf/filtered_665_accession_snp_maf0.01_gene.vcf.gz --chr 4 --recode --recode-INFO-all --stdout | bgzip -c > vcf/filtered_665_accession_snp_maf0.01_gene_Chr4.vcf.gz &
vcftools --gzvcf vcf/filtered_665_accession_snp_maf0.01_gene.vcf.gz --chr 5 --recode --recode-INFO-all --stdout | bgzip -c > vcf/filtered_665_accession_snp_maf0.01_gene_Chr5.vcf.gz &

# generate SNP matrix
vcftools --gzvcf vcf/filtered_665_accession_snp_maf0.01_gene_Chr1.vcf.gz --012 --out vcf/filtered_665_accession_snp_maf0.01_gene_Chr1 &
vcftools --gzvcf vcf/filtered_665_accession_snp_maf0.01_gene_Chr2.vcf.gz --012 --out vcf/filtered_665_accession_snp_maf0.01_gene_Chr2 &
vcftools --gzvcf vcf/filtered_665_accession_snp_maf0.01_gene_Chr3.vcf.gz --012 --out vcf/filtered_665_accession_snp_maf0.01_gene_Chr3 &
vcftools --gzvcf vcf/filtered_665_accession_snp_maf0.01_gene_Chr4.vcf.gz --012 --out vcf/filtered_665_accession_snp_maf0.01_gene_Chr4 &
vcftools --gzvcf vcf/filtered_665_accession_snp_maf0.01_gene_Chr5.vcf.gz --012 --out vcf/filtered_665_accession_snp_maf0.01_gene_Chr5 &

# POCRE screen
/home/xfu/R/3.2.4/bin/Rscript POCRE.R vcf/filtered_665_accession_snp_maf0.01_gene_Chr1.012
/home/xfu/R/3.2.4/bin/Rscript POCRE.R vcf/filtered_665_accession_snp_maf0.01_gene_Chr2.012
/home/xfu/R/3.2.4/bin/Rscript POCRE.R vcf/filtered_665_accession_snp_maf0.01_gene_Chr3.012
/home/xfu/R/3.2.4/bin/Rscript POCRE.R vcf/filtered_665_accession_snp_maf0.01_gene_Chr4.012
/home/xfu/R/3.2.4/bin/Rscript POCRE.R vcf/filtered_665_accession_snp_maf0.01_gene_Chr5.012

P=0.05

cat vcf/filtered_665_accession_snp_maf0.01_gene_Chr*.012.pos.select.$P > vcf/filtered_665_accession_snp_maf0.01_gene.pos.select.$P
vcftools --gzvcf vcf/filtered_665_accession_snp_maf0.01_gene.vcf.gz --positions vcf/filtered_665_accession_snp_maf0.01_gene.pos.select.$P --recode --recode-INFO-all --stdout | bgzip -c > vcf/filtered_665_accession_snp_maf0.01_gene.pos.select.$P.vcf.gz
bcftools index -t vcf/filtered_665_accession_snp_maf0.01_gene.pos.select.$P.vcf.gz
bcftools index --nrecords vcf/filtered_665_accession_snp_maf0.01_gene.pos.select.$P.vcf.gz # get total variant count

# SNP from Ensembl
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-39/vcf/arabidopsis_thaliana/arabidopsis_thaliana.vcf.gz -P vcf
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-39/vcf/arabidopsis_thaliana/arabidopsis_thaliana.vcf.gz.tbi -P vcf
bcftools index --nrecords vcf/arabidopsis_thaliana.vcf.gz
vcftools --gzvcf vcf/arabidopsis_thaliana.vcf.gz --remove-indels --recode --recode-INFO-all --stdout | bgzip -c > vcf/arabidopsis_thaliana_snp.vcf.gz
bcftools index -t vcf/arabidopsis_thaliana_snp.vcf.gz
bcftools index --nrecords vcf/arabidopsis_thaliana_snp.vcf.gz

# add Ensembl variation ID
bcftools annotate --annotations vcf/arabidopsis_thaliana_snp.vcf.gz --columns ID --output vcf/filtered_665_accession_snp_maf0.01_gene.pos.select.$P.ensembl.vcf.gz --output-type z vcf/filtered_665_accession_snp_maf0.01_gene.pos.select.$P.vcf.gz
zcat vcf/filtered_665_accession_snp_maf0.01_gene.pos.select.$P.ensembl.vcf.gz|grep -v '#'|awk '{print $1"\t"$2-1"\t"$2"\t"$3}' > vcf/filtered_665_accession_snp_maf0.01_gene.pos.select.$P.ensembl.bed
bedtools intersect -a tair10_gene.bed -b vcf/filtered_665_accession_snp_maf0.01_gene.pos.select.$P.ensembl.bed -wa -wb|cut -f4,10 > tair10_gene.SetID.$P


# PLINK (https://www.cog-genomics.org/plink/2.0)
wget -c https://www.cog-genomics.org/static/bin/plink180612/plink_linux_x86_64.zip -P src
mkdir plink

#zcat vcf/filtered_665_accession_snp_maf0.01_gene.pos.select.ensembl.vcf.gz|head -1000|bgzip -c > test.vcf.gz
#zcat test.vcf.gz|grep -v '#'|awk '{print $1"\t"$2-1"\t"$2"\t"$3}' > test.tmp
#bedtools intersect -a tair10_gene.bed -b test.tmp -wa -wb|cut -f4,10 > test.SetID

# run SKAT on cluster
qsub qsub.pbs

bedtools window -a tair10_gene.bed -b tair10_gene.bed -w 20000|cut -f4,10 > tair10_gene_local_20kb

