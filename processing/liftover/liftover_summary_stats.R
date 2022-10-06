#!/usr/bin/env Rscript

## start like this (example with convert_markername)
# Rscript --vanilla liftover_summary_stats.R \
# --summary_stats example_Kunkle2020_ADGC_AA_META_Model1_SummaryStats.withAlleleFreqs.txt \
# --from_build "hg19" \
# --to_build "hg38" \
# --chr "Chr" \
# --pos "Pos" \
# --convert_markername "MarkerName" \
# --out example.txt

## start like this (example without convert_markername)
# Rscript --vanilla liftover_summary_stats.R \
# --summary_stats example_finngen_R5_G6_ALZHEIMER_EXMORE.txt \
# --from_build "hg38" \
# --to_build "hg19" \
# --chr "#chrom" \
# --pos "pos" \
# --out example2.txt

library("optparse")
 
option_list = list(
  make_option(c("-s", "--summary_stats"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-f", "--from_build"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-t", "--to_build"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-m", "--convert_markername"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-c", "--chr"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-p", "--pos"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
   make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

require(data.table)
require(dplyr)
require(bigsnpr) # install.packages("bigsnpr")
data <- fread(opt$summary_stats,header=T)
data2 <- rename(data, chr = opt$chr) %>% rename(pos = opt$pos)
lifted <- snp_modifyBuild(data2, 'insert_liftOver_binary_path', from = opt$from_build, to = opt$to_build)

if (!is.null(opt$convert_markername)) {
    markernames <- lifted %>% pull(opt$convert_markername)
    chr_part <- gsub(":.*", "", markernames)
    alleles_part <- gsub(".*:[0-9]+", "", markernames)
    lifted[,opt$convert_markername] <- paste(chr_part, ":", lifted$pos, alleles_part, sep="")
}

# convert chr and pos back to original names
lifted <- lifted %>% filter(!is.na(pos))
names(lifted)[names(lifted)=="chr"] <- opt$chr
names(lifted)[names(lifted)=="pos"] <- opt$pos
write.table(lifted, file=opt$out, quote=FALSE,row.names=F,col.names=T,sep="\t")
