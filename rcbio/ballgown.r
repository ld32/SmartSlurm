#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (sample groups, such as 1,1,1,2,2,2).\n", call.=FALSE)
}

.libPaths("/n/shared_db/misc/rcbio/rlib/4.4.2")

#source("http://bioconductor.org/biocLite.R")
#biocLite("ballgown", suppressUpdates=TRUE)


library(ballgown);  

bg = ballgown(dataDir='stringtieOut', samplePattern='sample', meas='all'); 
save(bg, file='ballgown.rda'); 

group=unlist(strsplit(args[1], split=","));  
group; 
pData(bg) = data.frame(id=sampleNames(bg), group=group); 

whole_tx_table = texpr(bg, 'all')
stat_results = stattest(bg, feature='transcript', meas='FPKM', covariate='group'); 

trans=cbind(whole_tx_table, stat_results);  
write.table(trans[order(trans$qval,decreasing=F),], file='transcripts.txt', col.names=NA,quote=FALSE); 

#top 10 transcipt
top10=head(trans[order(trans$qval,decreasing=F),],10)

#source('/opt/rcbio-1.6/bin/geneFigure.r'); 
for (i in 1:10) {
    
    pdf(paste('gene_', i, '_',top10[i,9], '_all_samples_FPKM.pdf', sep=''));
    plotTranscripts(top10[i,9], bg, 
        samples=substr(colnames(top10)[10+2*(1:length(group))],6,100), 
        meas='FPKM', colorby='transcript')

    dev.off(); 
    
    #pdf(paste('transcript_',top10[i,9], '_compare_group_FPKM.pdf', sep=''));
    # !!!!!!!! this function only works for genes with multiple transcripts !!!!!!!!!!!
    #plotMeans(top10[i,9], bg, groupvar='group', meas='FPKM', colorby='transcript'); 
    #geneFigure(top10[i,9], bg, groupvar='group', meas='FPKM', colorby='transcript'); 
    #dev.off(); 
} 

stat_results = stattest(bg, feature='gene', meas='FPKM', covariate='group'); 
gene_expression = gexpr(bg)
genes=cbind(gene_expression, stat_results); 
write.table(genes[order(genes$qval,decreasing=F),], file='genes.txt', col.names=NA,quote=FALSE); 

cat("Done. Exiting from R script ...\n");
