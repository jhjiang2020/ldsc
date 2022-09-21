#!/usr/bin/env Rscript
# read arguments
library("optparse")
 
option_list = list(
  make_option(c("-f", "--sumstats"), type="character", default=NULL, 
              help="full summary statistics", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-p", "--pval"), type="character", default=NULL, 
              help="column name for p-values if filtering is required [optional]", metavar="character"),
  make_option(c("-n", "--snp_col"), type="character", default="rsID", 
              help="column name for rsids [default= %default]", metavar="character")

); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt$sumstats)){
  print_help(opt_parser)
  stop("Summary statistics file must be supplied (input file).n", call.=FALSE)
}



# load functions
library(dplyr)
source("/well/ocallaghan/users/qhs307/R/utils_R/GWAS_functions.R")

# MAKE SURE bcftools IS INSTALLED
bcftools="bcftools"

# load gwas sumstats
gwas.hg19 <- data.table::fread(opt$sumstats, header=T)
stopifnot(opt$snp_col %in% colnames(gwas.hg19))
rsids <- gwas.hg19 %>% pull(opt$snp_col)

## Get rsids with bcftools query
vcf = "/well/ocallaghan/users/qhs307/raw_data/dbSNP/dbSNP.hg38.vcf.gz"

gwas.hg38.map <- rsid2pos(bcftools = bcftools, SNP_rsid=rsids, vcf=vcf)

## Note this will only retain the first posid if duplicated
gwas.hg38 <- gwas.hg19
gwas.hg38$posID = gwas.hg38.map$posid[match(rsids, gwas.hg38.map$rsid)]

if(is.null(opt$pval)){
  data.table::fwrite(gwas.hg38[ , eval(opt$snp_col):=NULL], file=paste0(opt$out,".gz"),quote=F, sep="\t")
}else{
  gwas.hg38$P <- gwas.hg38 %>% pull(eval(opt$pval)) %>% as.numeric
  gwas.hg38 <-  gwas.hg38 %>% filter(P > 1e-323) ## minimal np.float64 is about 5e-324
  data.table::fwrite(gwas.hg38[ , c(eval(opt$snp_col),eval(opt$pval)):=NULL], file=paste0(opt$out,".gz"),quote=F, sep="\t")
}


