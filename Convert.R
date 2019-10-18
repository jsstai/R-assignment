#!/usr/bin/env Rscript
##########
# INTRO: #
##########
# FIELDS:
# [1]: URL or path to genotable file
# [2]: URL or path to snppos file
# [3]: String of non-chromosomal chromosome categories e.g. unknown, multiple, etc.,
#      separated by commas: 
#      ex.:-  unknown,multiple
# [4]: String of genotype groups, e.g. maize, teosinte, etc., with each group separated 
#      by double-commas and each genotype in group separated by a single-comma, 
#      ex.:-  ZMMIL,ZMMLR,ZMMMR,,ZMPBA,ZMPIL,ZMPJA
#
# EXAMPLE USAGE:
# args = c("https://raw.githubusercontent.com/EEOB-BioData/BCB546X-Fall2019/master/assignments/UNIX_Assignment/fang_et_al_genotypes.txt", "https://raw.githubusercontent.com/EEOB-BioData/BCB546X-Fall2019/master/assignments/UNIX_Assignment/snp_position.txt", "unknown,multiple", "ZMMIL,ZMMLR,ZMMMR,,ZMPBA,ZMPIL,ZMPJA")
#

# load commandline arguments into string vector args
args = commandArgs(trailingOnly=TRUE)

####################
# INPUT TEST PHASE #
####################
# check if there's at least one argument; if not, give error and direct to longer help message:
if (length(args)==0) {
  stop("At least one argument must be supplied; use option -help or -h for more details", call.=FALSE)

# check if help arg requested, display help if so.
} else if (args[1]=="-help" | args[1]=="-h") {
  stop("\nUsage: Convert.R PATHofSNPDataSet.txt PATHofSNPPositionsFile.txt ExcludeTerm1,ET2,etc. Group1Term1,G1T2,G1T3,,G2T1,G2T2,etc.\n\n Example: Convert.R https://raw.githubusercontent.com/EEOB-BioData/BCB546X-Fall2019/master/assignments/UNIX_Assignment/fang_et_al_genotypes.txt https://raw.githubusercontent.com/EEOB-BioData/BCB546X-Fall2019/master/assignments/UNIX_Assignment/snp_position.txt unknown,multiple ZMMIL,ZMMLR,ZMMMR,,ZMPBA,ZMPIL,ZMPJA", call.=TRUE)

# test if correct number of args exist; if not, give error and direct to longer help message:
} else if (length(args)!=4) {
  stop("Incorrect number of arguments supplied; use option -help or -h for more details", call.=FALSE)
}

##############
# PREP PHASE #
##############
# Check-for and/or load required packages
if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)

# Read files into variables
# Read arg one and record with and without heading as the genotable file:
genotable <- read_tsv(args[1], col_names = TRUE)
# Read arg two and record as SNP position file:
snppos <- read_tsv(args[2])
# Read arg three and split into genotype grouping e.g. maize vs. teosinte
excludes <- unlist(strsplit(args[3], ","))
# Read arg four and split into separated genotype grouping strings e.g. maize vs. teosinte
groupings <- unlist(strsplit(args[4], ",,"))

##################
# MAIN FUNCTIONS #
##################
# UNLOOPED SNP LIST WORK
# Generate list of unique chromosomal classificaitons from snppos:
chroms <- unique(snppos$Chromosome)

# Generate two lists of the chroms present in snppos based on whether
# they are present or absent in the arg-provided exclude list e.g. unknown, multiple:
chroms_true <- chroms[!chroms %in% excludes]
chroms_false <- chroms[chroms %in% excludes]

# OUTER LOOP BEGINNING: Set up subgenotables for each genotype grouping
for (grouping in groupings[])
{
# Split grouping string into a separated group of genotypes for co-analysis
# e.g. all maize genotypes as a group, all teosinte genotypes as a group:
  group <- c(unlist(strsplit(grouping[], ",")))
  
# Create subgenotable for entire phenotypic grouping from input-field 3:
  subgenotable <- filter(genotable, `Group` %in% group)
  
# Create transposed subgenotable:
  subgenotable_t <- t(subgenotable)
  
# Set Sample IDs as column names:
  colnames(subgenotable_t) <- as.character(unlist(subgenotable_t[1,]))
  subgenotable_t = subgenotable_t[-1, ]
  
# Set row names as first column:
  subgenotable_t <- rownames_to_column(as.data.frame(subgenotable_t), "SNP_ID")
  
# INNER LOOP 1: Generate files for each of the excluded false chromosome categories:
# e.g. unknown, multiple, etc.
  for (chrom in chroms_false[]) {
# Generate filename for chrom:
    chromname <- paste("SNPs", grouping[1], "chrom", chrom[1], basename(args[1]), sep=".")
    
# Generate SNP list for chrom:
    SNPs <- select(filter(snppos, `Chromosome` %in% chrom), `SNP_ID`, `Chromosome`, `Position`)
    
# Generate unsorted output table by filtering subgenotable_t by SNPs and appending result to SNPs
    assign(chromname, as.data.frame(append(SNPs, filter(subgenotable_t, `SNP_ID` %in% SNPs[[1]]))))
  }

# INNER LOOP 2: Generate files for each of the true chromosome IDs:
# e.g. 1, 2, 4, etc.
  for (chrom in chroms_true[]) {
# Generate filenames for chrom:
    chromname_incr <- paste("SNPs", grouping[1], "chrom", chrom[1], "incr", basename(args[1]), sep=".")
    chromname_decr <- paste("SNPs", grouping[1], "chrom", chrom[1], "decr", basename(args[1]), sep=".")
    
# Generate SNP list for chrom:
    SNPs <- select(filter(snppos, `Chromosome` %in% chrom), `SNP_ID`, `Chromosome`, `Position`)
    
# Generate preliminary output table by filtering subgenotable_t by SNPs and appending result to SNPs
    outfile_alpha <- as.data.frame(append(SNPs, select(filter(subgenotable_t, `SNP_ID` %in% SNPs[[1]]),-1)))
# Convert "?/?"s to "-/-"s
    outfile_beta <- mutate_all(outfile_alpha, ~gsub("[?]/[?]","-/-",.))
  
# Generate the required output table that is ordered by increasing position
    assign(chromname_incr, arrange(outfile_alpha, `Position`))
    
# Generate the required output table that is ordered by decreasing position
    assign(chromname_decr, arrange(outfile_beta, desc(`Position`)))
  }
}

##########################
# CLEANUP & EXPORT PHASE #
##########################
# Remove intermediates
rm(args, genotable, snppos, excludes, groupings, chroms, chroms_true, chroms_false, grouping, group, subgenotable, subgenotable_t, chrom, chromname, chromname_incr, chromname_decr, SNPs, outfile_alpha, outfile_beta)

# Write output files
for (output in ls())
{
  write.table(output, file = output, sep = "\t", col.names = TRUE)
}

# Remove output intermediate
rm(output)

# Thanks for reading!
# Did you know that botanically speaking, pumpkin and rice are both types of fruit?