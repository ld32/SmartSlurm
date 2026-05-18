#!/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  stop("Usage: EdgeR-v0.1.r [aligner_output_path, such as bwaOut, bowtieOut or starOut]", call.=FALSE)
}


files=list.files(args[1],pattern="read.count.txt", recursive=TRUE)

if(length(files)==0) {
   stop("No count file is found in directory: ", args[1], call.=FALSE)
}


group=as.numeric(substr(files,6,6))

.libPaths("/n/shared_db/misc/rcbio/rlib/4.4.2") 

library(edgeR)

y=readDGE(files, path=args[1], group=group,labels=substr(dirname(files),1, 1000),header=FALSE)

write.table(y$count, file="rawCount.txt", quote=FALSE, sep="\t")

cat('Done. Raw count is in file rawCount.txt\n')

