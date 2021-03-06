
require(goseq)
require(data.table)
require(reshape)
require(ggplot2)
require(scales)
require(stringr)

## This code is more of less tailered for this set of analysis downloaded from Gene Lab.
## Therefore, attention needs to be paid on a few key variables:
## 1. adjusted p-value
## 2. log fold-change
## 3. Gene ID (ENSEMBL,REFSEQ)
## Therefore, you have to DYDD and make sure that these are accounted for before executing the scripts.
## This script is just for your reference
data = fread("GLDS-117_array_differential_expression.csv")
data = as.data.frame(data)

## This dataset has quite a bit of genes missing ENSEMBL ID
## You may skip this step
data = data[!is.na(data$ENSEMBL),]
nonunique = as.data.frame(table(data$ENSEMBL))
names(nonunique)=c("gene","count")
nonunique = nonunique[ nonunique$count>1, ]
data = data[ !data$ENSEMBL %in% nonunique$gene, ]
rownames(data) = data$ENSEMBL

## Quality Control
## You may skip the entire section
data.mean = data[,grep("Group.Mean",colnames(data))]
data.logFC = data[,grep("Log2fc",colnames(data))]
data.padj = data[,grep("Adj.p",colnames(data))]

## count number of sig genes by category
statistics = data.frame(regulation=c(),count=c(),comparison=c())
for(i in 1:ncol(data.logFC))
{
	temp = data.frame(logFC=data.logFC[,i],padj=data.padj[,i])
	sig.temp = temp[ temp$padj<0.1, ]
	stats = as.data.frame(table(sig.temp$logFC>0))
	if(nrow(stats)==0){
		next
	}
	colnames(stats)=c("regulation","count")
	stats$comparison = c(colnames(data.logFC)[i])
	statistics = rbind(statistics,stats)
}
statistics$comparison = gsub("Log2fc_","",statistics$comparison)
statistics$comparison = gsub("\\& ","",statistics$comparison)
statistics$comparison = gsub("Proton Radiation","P",statistics$comparison)
statistics$comparison = gsub("Sham-irradiated","S",statistics$comparison)
statistics$comparison = gsub(" day","D",statistics$comparison)
df = as.data.frame(matrix(unlist(str_split(statistics$comparison,"v")),ncol=2,byrow=TRUE))
names(df)=c("G1","G2")
statistics = cbind(statistics,df)

regtext = c("Down-regulated","Up-regulated")
statistics$regulation = as.numeric(statistics$regulation)
statistics$regulation=regtext[statistics$regulation]

png("overall_DEgenes_statistics_FDR10_bar.png")
ggplot(statistics,aes(x=comparison,y=count))+geom_bar(aes(fill=regulation),position="stack",stat="identity")+xlab("")+ylab("# DE Genes (FDR<0.1)")+theme_bw()+theme(panel.grid=element_blank())+scale_fill_manual(values=c("royalblue","orange"))+coord_flip()
dev.off()

## Plot lattice scatterplot
colnames(data.mean)=gsub("Group.Mean_","",colnames(data.mean))
colnames(data.mean)=gsub("[\\(\\)\\&]","",colnames(data.mean))
panel.cor <- function(x, y){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- round(cor(x, y), digits=2)
    txt <- paste0("R = ", r)
    cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}
# Customize upper panel
upper.panel<-function(x, y){
  points(x,y, pch = 19)
}
# Create the plots
png("paired_scatterplot.png",width=3000,height=3000,res=300)
pairs(data.mean, lower.panel = panel.cor,upper.panel = upper.panel)
dev.off()
#####
## END OF QC SECTION


### STARTING THE GOSEQ analysis
getGeneLists <- function(pwf, goterms, genome, ids){
  gene2cat <- getgo(rownames(pwf), genome, ids)
  cat2gene <- split(rep(names(gene2cat), sapply(gene2cat, length)),
                    unlist(gene2cat, use.names = FALSE))
  out <- list()
  for(term in goterms){
    tmp <- pwf[cat2gene[[term]],]
    tmp <- rownames(tmp[tmp$DEgenes > 0, ])
    out[[term]] <- data$SYMBOL[match(tmp,data$ENSEMBL)]
  }
  out
}

goseqRUN <- function(delist,status,genome)
{
	rownames(delist)=delist$gid
	isSigGene = delist$sig
	if(status == "down"){
		isSigGene[ delist$logFC>0 ] = FALSE
	}else{
		isSigGene[ delist$logFC<0 ] = FALSE
	}

	genes = as.integer(isSigGene)
	names(genes)=rownames(delist)

	pwf <- nullp(genes,genome,"ensGene")
	GO.wall = goseq(pwf,genome,"ensGene",test.cats=c("GO:BP"))
	##topGO.wall = head(GO.wall,30)
	topGO.wall = GO.wall[ GO.wall$over_represented_pvalue<0.01, ]
	goList = getGeneLists(pwf,topGO.wall$category,genome,"ensGene")
	topGO.wall$EnsemblID <- sapply(topGO.wall$category, function(x) paste0(goList[[x]], collapse = ","))
	return(topGO.wall)
}

goseqRUNfull <- function(delist,status,genome)
{
	rownames(delist)=delist$gid
	isSigGene = delist$sig
	if(status == "down"){
		isSigGene[ delist$logFC>0 ] = FALSE
	}else{
		isSigGene[ delist$logFC<0 ] = FALSE
	}

	genes = as.integer(isSigGene)
	names(genes)=rownames(delist)

	pwf <- nullp(genes,genome,"ensGene")
	GO.wall = goseq(pwf,genome,"ensGene",test.cats=c("GO:BP"))
	##topGO.wall = head(GO.wall,30)
	##topGO.wall = GO.wall[ GO.wall$over_represented_pvalue<0.01, ]
	topGO.wall = GO.wall
	goList = getGeneLists(pwf,topGO.wall$category,genome,"ensGene")
	topGO.wall$EnsemblID <- sapply(topGO.wall$category, function(x) paste0(goList[[x]], collapse = ","))
	return(topGO.wall)
}

testdata = data.frame(gid=rownames(data),logFC = data$`Log2fc_(Proton Radiation & 26 day)v(Proton Radiation & 3 day)`, p = data$`Adj.p.value_(Proton Radiation & 26 day)v(Proton Radiation & 3 day)`)
testdata$sig = testdata$p<0.1
head(testdata)

testdata.pathway.up = goseqRUN(testdata,"up","mm9")
testdata.pathway.down = goseqRUN(testdata,"down","mm9")

fulldata.pathway.up = goseqRUN(testdata,"up","mm9")
fulldata.pathway.down = goseqRUN(testdata,"down","mm9")
write.table(fulldata.pathway.up,file="GOBP.analysis.up.txt",sep="\t",quote=FALSE,row.names=F)
write.table(fulldata.pathway.down,file="GOBP.analysis.down.txt",sep="\t",quote=FALSE,row.names=F)

plotGOBAR <- function(GO.wall,direction)
{
	GO.wall$logP = -1*log(GO.wall$over_represented_pvalue,10)
	top10 = head(GO.wall,10)
	p1<-ggplot(top10,aes(x=reorder(term,logP),y=logP))+geom_bar(stat="identity",fill="royalblue")+theme_bw()+theme(panel.grid=element_blank())+coord_flip()+geom_hline(yintercept=seq(0,floor(max(top10$logP)/10)*10,by=10),color="white")+ylab("-log10(P)")+xlab("")
	p1<-p1+ggtitle(direction)
	return(p1)
}

fulldata.pathway.up.barplot = plotGOBAR(fulldata.pathway.up,"GOBP UP")
png("upregulated_GOBP_barplot.png",width=1500,height=800,res=300)
print(fulldata.pathway.up.barplot)
dev.off()

fulldata.pathway.down.barplot = plotGOBAR(fulldata.pathway.down,"GOBP DOWN")
png("downregulated_GOBP_barplot.png",width=1500,height=800,res=300)
print(fulldata.pathway.down.barplot)
dev.off()

## This part is for plotting of significant ontology terms
## This has to be adapted for different datasets
normdata = fread("GLDS-117_array_normalized-annotated.txt.gz")
colnames(normdata)=gsub("Mmus.C57.6T.lvCMC.","",colnames(normdata))
normdata = as.data.frame(normdata)
colnames(normdata)[ ncol(normdata)] = "shamIR.3d.postshamIR.rep1"
rownames(normdata)=normdata$ID
normdata = normdata[,-1]

grouping = data.frame(sample=colnames(normdata))
condition = as.data.frame(matrix(unlist(str_split(grouping$sample,"\\.")),ncol=4,byrow=TRUE))
names(condition)=c("condition","time","treatment","rep")
grouping = cbind(grouping,condition)
grouping$id = paste(grouping$condition,grouping$time)
grouping$id  = factor(grouping$id,levels=c("shamIR 1d","shamIR 3d","HZE 1d","HZE 3d", "HZE 5d", "HZE 12d", "HZE 26d"))

z = (normdata-rowMeans(normdata))/apply(normdata,1,sd)
plotdata = melt(as.matrix(z))
names(plotdata) = c("refid","sample","z")
plotdata$SYMBOL = data$SYMBOL[match(plotdata$refid,data$REFSEQ)]
plotdata = plotdata[!is.na(plotdata$SYMBOL),]
plotdata$group = grouping$id[match(plotdata$sample,grouping$sample)]
plotdata$group = factor(plotdata$group,levels=c("shamIR 1d","shamIR 3d","HZE 1d","HZE 3d", "HZE 5d", "HZE 12d", "HZE 26d"))
head(plotdata)

##complex heatmap
require(ComplexHeatmap)
library(circlize)

sig.z = z[ rownames(z) %in% data$REFSEQ[ rowSums(data.padj<0.1)>0],]

ha = HeatmapAnnotation(bar = sample(letters[1:3], 10, replace = TRUE), col = list(bar = c("a" = "red", "b" = "green", "c" = "blue")))

col_fun = colorRamp2(c(0, 3, 20), c("white", "lightblue", "royalblue"))	
heatcol = colorRamp2(c(-2, 0, 2), c("blue", "white", "orange"))
ha = HeatmapAnnotation(condition=grouping$condition,time=as.numeric(as.character(gsub("d","",grouping$time))),col=list(time=col_fun,condition=c("shamIR"="blue","HZE"="red")))
png("DEgenes_heatmap.png",width=1800,height=2500,res=300)
Heatmap(sig.z, clustering_distance_rows = "pearson", col = heatcol, show_row_names = FALSE, clustering_distance_columns = "pearson", clustering_method_rows = "ward.D2", clustering_method_columns="ward.D2", top_annotation=ha)
dev.off()

## start to plot here
require(scales)
require(ggplot2)

plotGOHeatmap <- function(genelist,pathname)
{
	tempdata = plotdata[ plotdata$SYMBOL %in% unlist(str_split(genelist,",")), ]
	p1 = ggplot(tempdata,aes(x=sample,y=SYMBOL))+geom_tile(aes(fill=z))+theme_bw()
	p1 = p1 + theme(panel.grid=element_blank(),panel.spacing=unit(0,"lines"),axis.ticks=element_blank())+xlab("")+ylab("")
	p1 = p1 + theme(axis.text.x=element_blank()) + facet_grid(.~group,space="free",scale="free")
	p1 = p1 + scale_fill_gradientn(values=rescale(c(-1,0,1)),colors=c("blue","white","orange"))
	pathname = gsub("\\:","",pathname)
	pathname = gsub(" ","_",pathname)
	filename = paste(pathname,"_GOheat.png",sep="")
	png(filename,width=2500,height=1200,res=300)
	print(p1)
	dev.off()
}



for(i in 1:nrow(testdata.pathway.down))
{
	plotGOHeatmap(testdata.pathway.down$EnsemblID[i],paste("Down",testdata.pathway.down$category[i],testdata.pathway.down$term[i]))
}

for(i in 1:nrow(testdata.pathway.up))
{
	plotGOHeatmap(testdata.pathway.up$EnsemblID[i],paste("Up",testdata.pathway.up$category[i],testdata.pathway.up$term[i]))
}
 
 