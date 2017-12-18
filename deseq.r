library(DESeq)
Args<-commandArgs()
count_file=read.table(Args[6],header=TRUE,row.names=1)
label=read.table(Args[7],row.names=1)
outdir=Args[8]
outdir2=Args[9]
condition=label[names(count_file),]
cds=newCountDataSet(count_file,condition)
cds=estimateSizeFactors(cds)

count_norm=counts(cds, normalized=TRUE)
dir=paste(Args[8],'.norm', sep="")
write.table(count_norm,file=dir,row.names=T,col.names=T,quote=F,sep='\t')

count_mat_norm=as.matrix(count_norm)
log2_count_mat_norm=log2(count_mat_norm+1)
dimnames(log2_count_mat_norm)=list(row.names(count_file),colnames(count_file))
dir=paste(Args[8],".log.norm", sep="")
write.table(log2_count_mat_norm,file=dir,row.names=T,col.names=T,quote=F,sep='\t')

cds=estimateDispersions(cds,sharingMode = "maximum",method='pooled',fitType="local")



library(RColorBrewer)
a = brewer.pal(11,"RdYlBu")[1:11]
b = brewer.pal(11,"PuOr")[1:11]
c = brewer.pal(11,"RdGy")[1:11]
n=0
factor_level=levels(condition)
while ( n<length(factor_level)-1 ) {
     m=n+1
     while ( m<length(factor_level) ) {
         output_csv_dir=paste(outdir2,"/",factor_level[n+1],"_vs_",factor_level[m+1],".com",sep="")
         output_pdf_dir=paste(outdir2,"/",factor_level[n+1],"_vs_",factor_level[m+1],".pdf",sep="")
         x = nbinomTest(cds,factor_level[n+1],factor_level[m+1])
         write.csv(x,output_csv_dir,row.names=FALSE, quote=FALSE)
         pdf(output_pdf_dir)

         plot(log2(x$baseMean),x$log2FoldChange,xlim=c(1,10),ylim=c(-6,6),xlab="mean reads/peak (log2)",ylab=paste("log2(",factor_level[m+1],"/",factor_level[n+1],")",sep=""),pch=20,cex=0.5,col=ifelse(x$padj<0.001,c[2],ifelse(x$padj<0.01,a[2],ifelse(x$padj<0.05,b[3],c[8]))))
         legend("topright",title="FDR",inset=0.05,c("0-0.001","0.001-0.01","0.01-0.05","0.05-1"),col=c(c[2],a[2],b[3],c[8]),pch=c(15,15,15),cex=1,bty="n")
         abline(h=1,lwd=1,lty=2,col="black")
         abline(h=-1,lwd=1,lty=2,col="black")
         dev.off()
         m=m+1
         }
     n=n+1
}
