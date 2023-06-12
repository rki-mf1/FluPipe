#!/usr/bin/env Rscript
#winterk
#read the minimap2 files per sample and segment in and give 5 best (in case of Number of residue matches) fitting refs back


library(data.table)
#library(parallel)
setDTthreads(threads = 10)
###############

#input path output path
args = commandArgs(trailingOnly=TRUE)

sample_in <- args[1]
path_out <- args[2]


########################################################################
#function to read the data in by segment
########################
read_minimap_output <- function(sample, path){
    #extract sample names
    sampleName <- gsub (".minimap.txt", "", basename(sample))
    #heck if the samples are empty, if yes return NA if nonempty return first hit
    if(file.size(sample) > 0){
      minimap_output <- fread(file=sample, header=FALSE, select=c(6,10), sep= "\t", fill=TRUE)
      #V6 Target sequence name
      #V10 Number of residue matches
      #collapse duplicated references and sum up the matching residues
      Ref_freq_matches<- aggregate(minimap_output$V10 ~ minimap_output$V6, data=minimap_output, FUN=sum)
      #order by matching residues
      ordered_Ref_freq_matches <- Ref_freq_matches[order(Ref_freq_matches[,2], decreasing = TRUE),]
      #if else in case we don't get 5 hit Refs
      if(nrow(ordered_Ref_freq_matches) >= 5){
        best_ref <-ordered_Ref_freq_matches[1:5,1]
      }
      else{
        best_ref <-ordered_Ref_freq_matches[1:nrow(ordered_Ref_freq_matches),1]
      }
      write.table(best_ref,file=paste0(path, "/", sampleName, "_best_refs.txt"), quote=FALSE, row.names = FALSE, col.names = FALSE)
    }
    #empty input empty output
    else{
      best_ref <- NULL
      write.table(best_ref,file=paste0(path, "/", sampleName, "_best_refs.txt"), quote=FALSE, row.names = FALSE, col.names = FALSE)
    }
}
read_minimap_output(sample_in, path_out)