source("/mnt/projects/wlwtan/cardiac_epigenetics/pipeline/wilsonNGS/RNAseq/bin/FGSEA.R")

args = commandArgs(TRUE)

defile=args[1]
dename=args[2]

de = read.table(defile,header=TRUE)
de$fpkm = rowMeans(de[,8:ncol(de)])
de = de[ de$fpkm>=1, ]
ranks = de$logFC
names(ranks)=de$hgnc
head(ranks)




FGSEA("/mnt/projects/wlwtan/cardiac_epigenetics/pipeline/wilsonNGS/RNAseq/testdata/SHAMTAC/pathDB/allpath.all.v7.1.symbols.Rdata",ranks,paste("ALLPATH_",dename,sep=""))
FGSEA("/mnt/projects/wlwtan/cardiac_epigenetics/pipeline/wilsonNGS/RNAseq/testdata/SHAMTAC/pathDB/GO.c5.all.v7.1.symbols.Rdata",ranks,paste("GO_",dename,sep=""))
FGSEA("/mnt/projects/wlwtan/cardiac_epigenetics/pipeline/wilsonNGS/RNAseq/testdata/SHAMTAC/pathDB/MSigDB.all.v7.1.symbols.Rdata",ranks,paste("MSigDB_",dename,sep=""))
FGSEA("/mnt/projects/wlwtan/cardiac_epigenetics/pipeline/wilsonNGS/RNAseq/testdata/SHAMTAC/pathDB/regulatory.c3.all.v7.1.symbols.Rdata",ranks,paste("Regulatory_",dename,sep=""))
FGSEA("/mnt/projects/wlwtan/cardiac_epigenetics/pipeline/wilsonNGS/RNAseq/testdata/SHAMTAC/pathDB/c2.cp.v7.1.symbols.Rdata",ranks,paste("KEGG_",dename,sep=""))



