require(goseq)
require(ggplot2)

args = commandArgs(TRUE)

defile = args[1]
outputname=args[2]

uptxt = paste(args[2],"_upGOBP.txt",sep="")
upplot = paste(args[2],"_upGOBP.png",sep="")
downtxt = paste(args[2],"_downGOBP.txt",sep="")
downplot = paste(args[2],"_downGOBP.png",sep="")

args[1]


delist = read.table(args[1],sep="\t",header=TRUE)
delist$fpkm = rowMeans(delist[,8:ncol(delist)])
delist$sig = delist$FDR<0.05 & delist$logFC>=0.5 & delist$fpkm>=1
head(delist)

rownames(delist)=gsub("\\..*","",rownames(delist))
isSigGene = delist$sig
table(isSigGene)

genes = as.integer(isSigGene)
names(genes)=rownames(delist)

pwf <- nullp(genes,"mm9","ensGene")
GO.wall = goseq(pwf,"mm9","ensGene",test.cats=c("GO:BP"))
head(GO.wall)
write.table(GO.wall,file=uptxt,sep="\t",quote=FALSE,row.names=F)
dev.off()

GO.wall$logP = -1*log(GO.wall$over_represented_pvalue,10)
top10 = head(GO.wall,10)
top10

png(upplot,width=2000,height=1000,res=300)
ggplot(top10,aes(x=reorder(term,logP),y=logP))+geom_bar(stat="identity",fill="royalblue")+theme_bw()+theme(panel.grid=element_blank())+coord_flip()+geom_hline(yintercept=seq(0,floor(max(top10$logP)/10)*10,by=10),color="white")+ylab("-log10(P)")+xlab("")
dev.off()




delist$sig = delist$FDR<0.05 & delist$logFC<=-0.5 & delist$fpkm>=1
head(delist)

isSigGene = delist$sig
table(isSigGene)

genes = as.integer(isSigGene)
names(genes)=rownames(delist)

pwf <- nullp(genes,"mm9","ensGene")
GO.wall = goseq(pwf,"mm9","ensGene",test.cats=c("GO:BP"))
head(GO.wall)
write.table(GO.wall,file=downtxt,sep="\t",quote=FALSE,row.names=F)
dev.off()

GO.wall$logP = -1*log(GO.wall$over_represented_pvalue,10)
top10 = head(GO.wall,10)
top10

png(downplot,width=2000,height=1000,res=300)
ggplot(top10,aes(x=reorder(term,logP),y=logP))+geom_bar(stat="identity",fill="royalblue")+theme_bw()+theme(panel.grid=element_blank())+coord_flip()+geom_hline(yintercept=seq(0,floor(max(top10$logP)/10)*10,by=10),color="white")+ylab("-log10(P)")+xlab("")
dev.off()

