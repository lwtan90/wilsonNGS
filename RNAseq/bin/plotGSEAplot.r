require(fgsea)

args = commandArgs(TRUE)

rdata = args[1]
pathname = args[2]

database="/mnt/projects/wlwtan/cardiac_epigenetics/pipeline/wilsonNGS/RNAseq/testdata/SHAMTAC/pathDB/allpath.all.v7.1.symbols.Rdata"
if(pathname == "ALLPATH"){
	database="/mnt/projects/wlwtan/cardiac_epigenetics/pipeline/wilsonNGS/RNAseq/testdata/SHAMTAC/pathDB/allpath.all.v7.1.symbols.Rdata"
}else if(pathname == "GO"){
	database="/mnt/projects/wlwtan/cardiac_epigenetics/pipeline/wilsonNGS/RNAseq/testdata/SHAMTAC/pathDB/GO.c5.all.v7.1.symbols.Rdata"
}else if(pathname == "MSIGDB"){
	database="/mnt/projects/wlwtan/cardiac_epigenetics/pipeline/wilsonNGS/RNAseq/testdata/SHAMTAC/pathDB/MSigDB.all.v7.1.symbols.Rdata"
}else if(database == "TF"){
	database="/mnt/projects/wlwtan/cardiac_epigenetics/pipeline/wilsonNGS/RNAseq/testdata/SHAMTAC/pathDB/regulatory.c3.all.v7.1.symbols.Rdata"
}


database

load(as.character(database))
load(rdata)

pathname = fgseaRes$pathway[ fgseaRes$padj<0.1 ]
pathname

for(i in pathname)
{
	plotfile = paste(rdata,i,"GSEAplot.pdf",sep="_")
	pdf(plotfile)
	print(plotEnrichment(path[[ as.character(i) ]],ranks))
	dev.off()

}
