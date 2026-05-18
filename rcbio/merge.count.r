#!/bin/env Rscript

.libPaths("/n/shared_db/misc/rcbio/rlib/4.4.2")

#source("http://bioconductor.org/biocLite.R")
#biocLite("edgeR", suppressUpdates=TRUE)

args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  stop("Usage: merge.count.r [aligner_output_path, such as bwaOut, bowtieOut or starOut]", call.=FALSE)
}


files=list.files(args[1],pattern="read.count.txt", recursive=TRUE)

if(length(files)==0) {
   stop("No count file is found in directory: ", args[1], call.=FALSE)
}

group=as.numeric(substr(files,6,6))

library(edgeR)
y=readDGE(files, path=args[1], group=group,labels=files,header=FALSE)
write.table(y$count, file="rawCount.txt", col.names=NA, quote=FALSE);

cat('Done. Raw count is in file rawCount.txt\n')
