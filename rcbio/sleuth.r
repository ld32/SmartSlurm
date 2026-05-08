#!/usr/bin/env Rscript

args<-commandArgs(TRUE)

sp=args[1]

#source("http://bioconductor.org/biocLite.R")

#biocLite("devtools", suppressUpdates=TRUE)

#biocLite("pachterlab/sleuth", suppressUpdates=TRUE)

#biocLite("biomaRt", suppressUpdates=TRUE)

.libPaths("/n/shared_db/misc/rcbio/rlib/4.4.2") 

library("sleuth")

sample=read.table("sample.lst", header = TRUE, stringsAsFactors=FALSE)

if (sp != "") {
    library("biomaRt")   
    mart=biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset=sp, host = "https://useast.ensembl.org")
    t2g=biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart)

    t2g=dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

    sout=sleuth_prep(sample, ~ group, target_mapping = t2g)
} else {
    sout=sleuth_prep(sample, ~ group)
}
    
sout=sleuth_fit(sout)

sout=sleuth_fit(sout, ~1, 'reduced')

sout=sleuth_lrt(sout, 'reduced', 'full')

# see above another test you can use is sleuth_wt. Here is the discussion on which one you should choose
# https://gist.github.com/jaquol/03f41f57dc6b0eacef101e9920f24d78

results=sleuth_results(sout, 'reduced:full', test_type = 'lrt')

results_ordered=results[order(results$qval),]

write.table(results_ordered, file='sleuth.DE_transcripts.txt', sep="\t",row.names=F, quote=F)

write.table( subset(results_ordered, qval <= 0.05), file='sleuth.DE_transcripts.qval_0.05.txt', sep="\t",row.names=F, quote=F)
