
require(fgsea)

FGSEA <- function(db,ranks,dbtype){
##depending on where you are working with mouse or human
##load("pathDB/Mm.c5.symbol.Rdata")
##full path to DB
## the variable stored is path
## we use hgnc remember

pathRdata = paste(dbtype,".RData",sep="")
pathTxt = paste(dbtype,"_FGSEA.txt",sep="")


load(db)
ls()

##data = read.table("de_SHAM_TAC.txt",header=TRUE)
##dim(data)


## This is to test for enrichment
##ranks = data$logFC
##names(ranks)=data$hgnc
pathgene=c()

fgseaRes <- fgsea(path,ranks,minSize=15,maxSize=500,nperm=1000)
sigPATH = fgseaRes[ order(padj, -abs(NES)), ]
for(i in 1:nrow(sigPATH)){ temp=paste0(unlist(sigPATH[i,8]),collapse=","); pathgene=c(pathgene,temp) }

sigPATH = sigPATH[,1:7]
sigPATH$gene = pathgene
print(head(sigPATH))

#sigPATH = as.data.frame(sigPATH[ sigPATH$padj<0.05, ])


##head(sigPATH)

## If you do GO, this will be GO

save(fgseaRes,file=pathRdata)
write.table(sigPATH,file=pathTxt,sep="\t",quote=FALSE,row.names=F)

}
