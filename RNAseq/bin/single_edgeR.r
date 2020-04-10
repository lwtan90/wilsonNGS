require(RUVSeq)
require(ggplot2)
require(edgeR)
require(EDASeq)
require(RColorBrewer)


args=commandArgs(TRUE)
group1 = args[1]
group2 = args[2]

outputfile = paste("de_",group1,"_",group2,".txt",sep="")
fpkmfile = paste("fpkm_",group1,"_",group2,".txt",sep="")
normfile = paste("norm_",group1,"_",group2,".txt",sep="")

bedup = paste("de_",group1,"_",group2,"_targetup.bed",sep="")
beddown = paste("de_",group1,"_",group2,"_targetdown.bed",sep="")
bedall = paste("de_",group1,"_",group2,"_all.bed",sep="")

genelength = read.table("/mnt/projects/wlwtan/cardiac_epigenetics/pipeline/rnaseq_pipeline/resource/mm9/mm9_genelength.txt")
names(genelength)=c("id","length")
gene2biotype = read.table("/mnt/projects/wlwtan/cardiac_epigenetics/pipeline/rnaseq_pipeline/resource/mm9/mm9_gene2biotype.txt")
names(gene2biotype)=c("id","hgnc","biotype")


filelist = read.table("countfile.txt")
names(filelist)=c("filename","sample","group")

filelist = filelist[ filelist$group %in% c(group1,group2), ]
filelist$group = factor(filelist$group,levels=c(group1,group2))
filelist$color = as.character(rep("blue"))
filelist$color[ filelist$group==group2 ] = "red"


## storing the count data as a DGEList object  
data = readDGE(file=filelist$filename,group=filelist$group,labels=filelist$sample)

## Filter away genes with very low counts across the samples.
## Assumption: basal transcription could be noise.
## Genes are dropped if they can't possibily be expressed in all samples. The count of these genes are 0 throughout.
## I followed the recommended number of 5-10 count per gene as a threshold in each library. 

filter = apply(data$counts,1,function(x) length(x[x>5])>=2)
filtered = data$counts[ filter, ]

## not sure if this is necessary as this will remove those no feature counts reported by htseq-count
filtered = filtered[grep("^_",rownames(filtered),invert=TRUE),]
set = newSeqExpressionSet(as.matrix(filtered),phenoData=data.frame(filelist$group,row.names=filelist$sample))
set = betweenLaneNormalization(set,which="upper")
countdata = normCounts(set)
fpkm = rpkm(countdata,gene.length=genelength$length)


## Normalization
## 1. sequencing depth
## 2. rna composition
## TMM normalization
design = model.matrix(~0+filelist$group)
design

y = DGEList(counts = as.matrix(filtered), group=filelist$group )
y = calcNormFactors(y)
y = estimateDisp(y, design)

#set = betweenLaneNormalization(set,which="upper")

write.table(filtered,file="raw_count.txt",sep="\t",quote=FALSE)
write.table(fpkm,file=fpkmfile,sep="\t",quote=FALSE)
write.table(countdata,file=normfile,sep="\t",quote=FALSE)


##design = model.matrix(~0+filelist$group,data=pData(set))
##y = DGEList(counts=counts(set),group=filelist$group)
##y = calcNormFactors(y,method="upperquartile")
##y = estimateDisp(y, design)
##y = estimateGLMCommonDisp(y,design)
##y = estimateGLMTagwiseDisp(y,design)


####This part should be editted according to the nature of the samples / experimental design
fit = glmFit(y,design)


##This step is to perform differential gene expression analysis
# I am using loglikihood test instead
lrt21 = glmLRT(fit,contrast=c(-1,1))
k21 = as.data.frame(topTags(lrt21,n=10000000))
k21$hgnc = gene2biotype$hgnc[match(rownames(k21),gene2biotype$id)]
k21$biotype = gene2biotype$biotype[match(rownames(k21),gene2biotype$id)]
k21 = cbind(k21,fpkm[match(rownames(k21),rownames(fpkm)),])
write.table(k21,file=outputfile,sep="\t",quote=FALSE)


##This step is to perform motif analysis
delist = k21
delist$fpkm = rowMeans(delist[,8:ncol(delist)])
delist = delist[ delist$FDR<0.05 & abs(delist$logFC)>=1 & delist$fpkm>=1,]

down.target = delist$hgnc[ delist$logFC<0 ]
up.target   = delist$hgnc[ delist$logFC>0 ]


gene = read.table("/mnt/projects/wlwtan/cardiac_epigenetics/pipeline/rnaseq_pipeline/resource/mm9/mm9_longest_tss.bed")
target.gene = gene[ gene$V4 %in% down.target, 1:4]
write.table(target.gene,file=beddown,sep="\t",quote=FALSE,row.names=F,col.names=F)
down = target.gene

target.gene = gene[ gene$V4 %in% up.target, 1:4]
write.table(target.gene,file=bedup,sep="\t",quote=FALSE,row.names=F,col.names=F)
up = target.gene

allgene = gene[ gene$V4 %in% delist$hgnc[ delist$fpkm>=1 ], ]
write.table(allgene,file=bedall,sep="\t",quote=FALSE,row.names=F,col.names=F)

system(paste("mkdir upmotif_",group1,group2,";   cd upmotif_",group1,group2,";  findMotifsGenome.pl ../",bedup,"   mm9 . -bg ../",bedall," -nomotif",sep=""),intern=TRUE)
system(paste("mkdir downmotif_",group1,group2,";   cd downmotif_",group1,group2,";  findMotifsGenome.pl ../",beddown,"   mm9 . -bg ../",bedall," -nomotif",sep=""),intern=TRUE)



require(ggplot2)

pcaplot_global = paste("globalpca_",group1,"_",group2,".png",sep="")
clusterplot_global = paste("globalclustering_",group1,"_",group2,".png",sep="")
pcaplot_degene = paste("degenepca_",group1,"_",group2,".png",sep="")
heatmapdegene = paste("degene_heatmap_",group1,"_",group2,".png",sep="")
scatterplot_notext = paste("scatterplot_notext_",group1,"_",group2,".png",sep="")
volcanoplot_notext = paste("volcanoplot_notext_",group1,"_",group2,".png",sep="")
scatterplot_withtext = paste("scatterplot_withtext_",group1,"_",group2,".png",sep="")
volcanoplot_withtext = paste("volcanoplot_withtext_",group1,"_",group2,".png",sep="")




control = fpkm[,colnames(fpkm) %in% filelist$sample[ filelist$group==group1] ]
condition2 = fpkm[,colnames(fpkm) %in% filelist$sample[ filelist$group==group2] ]

statistics = data.frame(control=rowMeans(control),condition2=rowMeans(condition2))

k21$fpkm = rowMeans(k21[,8:ncol(k21)])

k21$sig = k21$FDR<0.05 & abs(k21$logFC)>log(1.5,2) & k21$fpkm>=1
k21$up = k21$FDR<0.05 & k21$logFC>log(1.5,2) & k21$fpkm>=1

statistics$k21 = rownames(statistics)%in%rownames(k21)[ k21$sig ]
statistics$k21.up = rownames(statistics)%in%rownames(k21)[ k21$up ]
statistics$hgnc = k21$hgnc[match(rownames(statistics),rownames(k21))]



###For full data
fpkm = log(fpkm+1,2)
fpkm = (fpkm-rowMeans(fpkm))/apply(fpkm,1,sd)

pca = prcomp(t(fpkm),scale=TRUE,center=TRUE)
eigs <- pca$sdev^2
percentage_explained = eigs/sum(eigs)*100
percentage_explained

pcadata = as.data.frame(pca$x)
pcadata$group = filelist$group[match(rownames(pcadata),filelist$sample)]

pcaxl = paste("PC1(",format(round(percentage_explained[1],2),2),"%)",sep="")
pcayl = paste("PC2(",format(round(percentage_explained[2],2),2),"%)",sep="")

png(pcaplot_global,width=1500,height=2000,res=300)
ggplot(pcadata,aes(x=PC1,y=PC2))+geom_point(aes(color=group),shape=18,size=1)+theme_bw()+theme(panel.grid=element_blank())+xlab(pcaxl)+ylab(pcayl)+scale_color_manual(values=c("blue","red"))
dev.off()


hc = hclust(as.dist(1-cor(fpkm,method="spearman")),method="ward.D2")
png(clusterplot_global,width=2500,height=2000,res=300)
plot(hc)
dev.off()

####Only for DE genes
fpkm.sig = fpkm[ rownames(fpkm) %in% rownames(k21)[ k21$sig ], ]
fpkm.sig = (fpkm.sig-rowMeans(fpkm.sig))/apply(fpkm.sig,1,sd)

pca = prcomp(t(fpkm.sig),scale=TRUE,center=TRUE)
eigs <- pca$sdev^2
percentage_explained = eigs/sum(eigs)*100
percentage_explained

pcadata = as.data.frame(pca$x)
pcadata$group = filelist$group[match(rownames(pcadata),filelist$sample)]

pcaxl = paste("PC1(",format(round(percentage_explained[1],2),2),"%)",sep="")
pcayl = paste("PC2(",format(round(percentage_explained[2],2),2),"%)",sep="")

png(pcaplot_degene,width=1500,height=2000,res=300)
ggplot(pcadata,aes(x=PC1,y=PC2))+geom_point(aes(color=group),shape=18,size=1)+theme_bw()+theme(panel.grid=element_blank())+xlab(pcaxl)+ylab(pcayl)+scale_color_manual(values=c("blue","orange"))
dev.off()

filelist = filelist[match(colnames(fpkm),filelist$sample),]

require(gplots)

hc = hclust(as.dist(1-cor(fpkm.sig,method="spearman")),method="ward.D2")
hr = hclust(as.dist(1-cor(t(fpkm.sig),method="pearson")),method="ward.D2")

fpkm.sig = as.matrix(fpkm.sig)
rownames(fpkm.sig) = k21$hgnc[match(rownames(fpkm.sig),rownames(k21))]

png(heatmapdegene,width=1500,height=2300,res=300)
heatmap.2(as.matrix(fpkm.sig),Colv=as.dendrogram(hc),Rowv=as.dendrogram(hr),scale="none",breaks=c(seq(-2,2,by=0.04)),col=colorRampPalette(c("blue","white","orange"))(100),density.info="none",trace="none",ColSideColors=as.character(filelist$color))
dev.off()


sigcondition2 = statistics[ statistics$k21, ]
table(sigcondition2$k21.up)
head(sigcondition2)
head(statistics)
png(scatterplot_notext,width=2000,height=2000,res=300)
ggplot(statistics,aes(x=control+1,y=condition2+1))+geom_point(color="grey60",shape=18,size=1)+geom_point(aes(x=control+1,y=condition2+1),color="red",size=1.25,shape=18,data=sigcondition2)+theme_bw()+theme(panel.grid=element_blank())+xlab(paste("FPKM of ",group1))+ylab(paste("FPKM of ",group2))+scale_x_log10()+scale_y_log10()
dev.off()

png(scatterplot_withtext,width=2000,height=2000,res=300)
ggplot(statistics,aes(x=control+1,y=condition2+1))+geom_point(color="grey60",shape=18,size=1)+geom_point(aes(x=control+1,y=condition2+1),color="red",size=1.25,shape=18,data=sigcondition2)+geom_text(aes(x=control+1,y=condition2+1,label=hgnc),data=sigcondition2)+theme_bw()+theme(panel.grid=element_blank())+xlab(paste("FPKM of ",group1))+ylab(paste("FPKM of ",group2))+scale_x_log10()+scale_y_log10()
dev.off()


k21$logP = log(k21$PValue,10)*-1
de = k21[ k21$sig, ]

png(volcanoplot_notext,width=2000,height=2000,res=300)
ggplot(k21,aes(x=logFC,y=logP))+geom_point(aes(color=sig),shape=18,size=1)+theme_bw()+theme(panel.grid=element_blank())+xlab(paste("log2FC(",group2,"/",group1,")",sep=""))+ylab("-log10(P-Value)")+scale_color_manual(values=c("grey60","red"))
dev.off()

png(volcanoplot_withtext,width=2000,height=2000,res=300)
ggplot(k21,aes(x=logFC,y=logP))+geom_point(aes(color=sig),shape=18,size=1)+geom_text(aes(x=logFC+0.2,y=logP,label=hgnc),data=de)+theme_bw()+theme(panel.grid=element_blank())+xlab(paste("log2FC(",group2,"/",group1,")",sep=""))+ylab("-log10(P-Value)")+scale_color_manual(values=c("grey60","red"))
dev.off()




