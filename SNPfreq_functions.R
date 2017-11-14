# get minor allele frequency
minorAlleleFrequency <- function(SNPmat){
  # get minor allele frequency 
  # from matrix when SNPs are in columns
  maf <- colMeans(SNPmat, na.rm = TRUE)/2
  # find maf > 0.5 (common alleles)
  common <- which(maf > 0.5)
  # fix maf
  maf[common] <- 1 - maf[common]
  return(maf)
}

# subset the dataset to only keep SNPs with MAF > 0.05
filterSNPs <- function(SNPmat, minMAF = 0.05){
  # function to filter a SNP matrix
  # to only keep SNPs with maf > 0.05
  maf <- minorAlleleFrequency(SNPmat)
  SNPmat <- SNPmat[,maf > minMAF]
  return(SNPmat)
}

# function to get hist of minor allele freq
mafHist <- function(SNPmat){
  maf <- minorAlleleFrequency(SNPmat)
  hist(maf)
}
