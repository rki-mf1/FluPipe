#!/usr/bin/env Rscript
#winterk
#get variants with good qualities that failed the variant filtering due to the lofreq strand bias filter


library(data.table)
library(R.utils, warn.conflicts = FALSE) #.gz input
setDTthreads(threads = 10)
###############

#input path output path
args = commandArgs(trailingOnly=TRUE)

file_in <- args[1]
path_out <- args[2]

filter_special_variant_fail_case <- function(sample_name, path){
  #get sample name for output file name
  sampleName <- gsub (".no.filter.filtered.vcf.gz", "", basename(sample_name))
  #get the variant list per sample
  vcf <- fread(file=sample_name, header=TRUE, skip=14)
  #only keept good quality variants
  goodqual <- vcf[which(vcf$QUAL>=10000) ,]
  #only keep the failed ones 
  goodqual_notPASSED <- goodqual[which(goodqual$FILTER != "PASS"),]
  
  #in case there are no variants - done, write output
  if(nrow(goodqual_notPASSED) == 0){
    bed_tmp <- NULL
    write.table(bed_tmp, file=paste0(path, "/", sampleName,".special_case_variant_mask.bed"), quote = FALSE,  sep ="\t", col.names = FALSE, row.names = FALSE )
    
  }
  else{
    #more filtering
    #split INFO field to get the allele freq value
    goodqual_notPASSED[, c("DP", "AF", "SB", "DP4") := tstrsplit(INFO, ";", fixed=TRUE)]
    goodqual_notPASSED$AF2<-gsub("AF=","",as.character(goodqual_notPASSED$AF))
    goodqual_notPASSED$AF2<-as.numeric(goodqual_notPASSED$AF2)
    #filter variants with allele frequency >= 0.9
    goodqual_notPASSED_AFhigh <-goodqual_notPASSED[which(goodqual_notPASSED$AF2 >= 0.9),]
    
    #prepare data for the masking bed file
    #variant position minus one because in bed files coordinates start with 0
    goodqual_notPASSED_AFhigh$start <- as.data.frame(goodqual_notPASSED_AFhigh[,2])-1
    #get name, start and stop position
    bed_tmp <- cbind(goodqual_notPASSED_AFhigh[,1], goodqual_notPASSED_AFhigh$start, goodqual_notPASSED_AFhigh[,2])
    
    
    write.table(bed_tmp, file=paste0(path, "/", sampleName, ".special_case_variant_mask.bed"), quote = FALSE,  sep ="\t", col.names = FALSE, row.names = FALSE )
  }
  
}

filter_special_variant_fail_case(file_in, path_out)