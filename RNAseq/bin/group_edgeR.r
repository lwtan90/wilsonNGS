require(edgeR)
require(RUVSeq)

filelist = read.table("countfile.txt")
names(filelist)=c("countfile","sample","group")

data = readDGE(file=filelist$countfile,group=filelist$group,labels=filelist$sample)
filter = apply(data$counts,1,function(x) length(x[x>5])>=2)
filtered = data$counts[ filter, ]

metatag = grep("^__",rownames(filtered))
filtered = filtered[-metatag,]
summary(filtered)

y = DGEList(counts=filtered,group=filelist$group)
y = calcNormFactors(y)
group = filelist$group
group = factor(group,levels=c("NC_W8","NC_W16","HFHS_W8","HFHS_W16"))
group

design = model.matrix(~0+group)
design

y = estimateDisp(y,design)
fit = glmQLFit(y,design)

genelength = read.table("/mnt/projects/wlwtan/cardiac_epigenetics/pipeline/rnaseq_pipeline/resource/mm9/mm9_genelength.txt")
names(genelength)=c("id","length")

set = newSeqExpressionSet(as.matrix(filtered),phenoData=data.frame(filelist$group,row.names=filelist$sample))
set = betweenLaneNormalization(set,which="upper")
countdata = normCounts(set)
fpkm = rpkm(countdata,gene.length=genelength$length)
write.table(fpkm,file="fpkm.txt",sep="\t",quote=FALSE)


gene2biotype = read.table("/mnt/projects/wlwtan/cardiac_epigenetics/pipeline/rnaseq_pipeline/resource/mm9/mm9_gene2biotype.txt")
names(gene2biotype)=c("id","hgnc","biotype")

qlf = as.data.frame(topTags(glmQLFTest(fit,contrast=c(-1,1,0,0)),n=100000000))
qlf$hgnc = gene2biotype$hgnc[match(rownames(qlf),gene2biotype$id)]
qlf$biotype = gene2biotype$biotype[match(rownames(qlf),gene2biotype$id)]
qlf = cbind(qlf,fpkm[match(rownames(qlf),rownames(fpkm)),])
summary(qlf)
write.table(qlf,file="de_NCw8_NCw16.txt",sep="\t",quote=FALSE)

qlf = as.data.frame(topTags(glmQLFTest(fit,contrast=c(-1,0,1,0)),n=100000000))
qlf$hgnc = gene2biotype$hgnc[match(rownames(qlf),gene2biotype$id)]
qlf$biotype = gene2biotype$biotype[match(rownames(qlf),gene2biotype$id)]
qlf = cbind(qlf,fpkm[match(rownames(qlf),rownames(fpkm)),])
summary(qlf)
write.table(qlf,file="de_NCw8_HFHSw8.txt",sep="\t",quote=FALSE)

qlf = as.data.frame(topTags(glmQLFTest(fit,contrast=c(0,-1,0,1)),n=100000000))
qlf$hgnc = gene2biotype$hgnc[match(rownames(qlf),gene2biotype$id)]
qlf$biotype = gene2biotype$biotype[match(rownames(qlf),gene2biotype$id)]
qlf = cbind(qlf,fpkm[match(rownames(qlf),rownames(fpkm)),])
summary(qlf)
write.table(qlf,file="de_NCw16_HFHSw16.txt",sep="\t",quote=FALSE)

qlf = as.data.frame(topTags(glmQLFTest(fit,contrast=c(1,-1,-1,1)),n=100000000))
qlf$hgnc = gene2biotype$hgnc[match(rownames(qlf),gene2biotype$id)]
qlf$biotype = gene2biotype$biotype[match(rownames(qlf),gene2biotype$id)]
qlf = cbind(qlf,fpkm[match(rownames(qlf),rownames(fpkm)),])
summary(qlf)
write.table(qlf,file="de_HFHSw8_HFHSw16.txt",sep="\t",quote=FALSE)




