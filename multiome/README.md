# Wilson's guide to 10x Multiome GEX + ATAC analysis  
Date: 12/2/2021  
Author: Wilson  

## Description  
I use R package Seurat to run most of the analysis here.
This tutorial assumes that you have generated the output from cellranger arc.
If not, I will create additional tutorial to morph other possible input for the code.
Right now, the tutorial assumes the input files per sample to be:  
- matrix h5 (filtered_feature_bc_matrix.h5)  
- fragment file (fragments.tsv.gz)  


## Project 1: Differential GEX + ATAC analysis between two conditions  
Loading required packages:  

```
library(Signac) ## needed for single-cell chromatin analysis
library(Seurat) ## base package for single-cell analysis 
library(ggplot2) ## for plotting
library(JASPAR2020) ## needed for ChromVAR analysis
library(TFBSTools) ## needed for ChromVAR analysis
library(EnsDb.Rnorvegicus.v79) ## needed for one function (dont have to be exact, unless you are running chromvar)
library(BSgenome.Rnorvegicus.UCSC.rn6) ## needed for one function (dont have to be exacr, unless you are running chromvar)
```  

### Preparing annotation  
In the original tutorial by Sajita's lab, the loading of annotation seems straightforward. You just require the standard package from EnsDb and BSgenome for the organism of interest. I have listed the relevant packages below (you can find the rest of the organism online):  
| Organism | Genome Build | EnsDb | BSGenome | Tested by Myself |
| --- | --- | --- | --- | --- |
| Human | hg38 | EnsDb.Hsapiens.v86 | BSgenome.Hsapiens.UCSC.hg38 | No |
| Human | hg19 | EnsDb.Hsapiens.v75 | BSgenome.Hsapiens.UCSC.hg19 | Yes |
| Mouse | mm10 | EnsDb.Mmusculus.v79 | BSgenome.Mmusculus.UCSC.mm10 | Yes |
| Mouse | mm9 | EnsDb.Mmusculus.v75 | BSgenome.Mmusculus.UCSC.mm9 | Yes |
| Rat | rn6 | NA | BSgenome.Rnorvegicus.UCSC.rn6 | Yes |
| Rat | rn5 | EnsDb.Rnorvegicus.v79 | BSgenome.Rnorvegicus.UCSC.rn5 | Yes |  

### Reading 10x HDF5 input  
This function takes in h5 file, fragment file, and sample label. It returns Seurat object with RNA and ATAC components.
This function also assume that you have prepped the annotation variable in the step earlier.
Note: Only run this function if you have multiome (not multi-omics, meaning the ATAC and GEX are generated separately).  
```
readARC <- function(h5file, fragfile, sample)
{
        # load the RNA and ATAC data
        #"filtered_feature_bc_matrix.h5"
        #"atac_fragments.tsv.gz"
        counts <- Read10X_h5(h5file)
        fragpath <- fragfile
        obj <- CreateSeuratObject(
                counts = counts$`Gene Expression`,
                assay = "RNA"
        )
        # create ATAC assay and add it to the object
        obj[["ATAC"]] <- CreateChromatinAssay(
                counts = counts$Peaks,
                sep = c(":", "-"),
                fragments = fragpath,
                annotation = annotation
        )
        DefaultAssay(obj) = "ATAC"
        obj <- NucleosomeSignal(obj)
        obj <- TSSEnrichment(obj)

	## What are each of these
        png(paste(sample,"QC_violin_seurat.png",sep="_"),width=3000,height=2000,res=300)
        p1 <- VlnPlot(object=obj,features=c("nCount_RNA", "nFeature_RNA", "nCount_ATAC", "nFeature_ATAC","TSS.enrichment", "nucleosome_signal"),ncol=3,pt.size=0)
        print(p1)
        dev.off()

	## The numbers here were set arbitrarily
	## Better check the violin plots prior to setting numbers
	## Questions: Does it matter if you set different threshold for different samples?
        obj <- subset(
                x = obj,
                subset = nCount_ATAC < 30000 &
                nCount_RNA < 20000 &
                nCount_ATAC > 1000 &
                nCount_RNA > 1000 &
                nucleosome_signal < 2 &
                TSS.enrichment > 1
        )

	## Perform peak-calling using macs after removing noise
        peaks <- CallPeaks(obj, macs2.path = "macs2")

        # remove peaks on nonstandard chromosomes and in genomic blacklist regions
        peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")

        # quantify counts in each peak
        macs2_counts <- FeatureMatrix(
                fragments = Fragments(obj),
                features = peaks,
                cells = colnames(obj)
        )

        # create a new assay using the MACS2 peak set and add it to the Seurat object
        obj[["peaks"]] <- CreateChromatinAssay(
                counts = macs2_counts,
                fragments = fragpath,
                annotation = annotation
        )

        DefaultAssay(obj) = "RNA"
        obj$samplename = rep(sample)
        return(obj)
}

```  


