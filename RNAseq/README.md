## RNA-seq pipelines  
<br/>  
<p align="center">
  <img height="600" src="https://github.com/lwtan90/wilsonNGS/blob/master/RNAseq/img/RNAseq_workflow.jpg">
</p>  
 
<br/>
  
This pipeline will take in paired-end FASTQ files, and produce gene count as output.  
This is a skeleton for basic analysis.  
There are various added analysis described at the bottom of the page.  
See [bin/rnaseq.sh](https://github.com/lwtan90/wilsonNGS/blob/master/RNAseq/bin/rnaseq.sh)  


```  

#!/bin/bash

##INPUT ARGUMENT (DO NOT EDIT)
##BASH DEFAULT USER ARGUMENT READING FROM COMMAND LINE
SAMPLENAME=$1
FASTQ_READ1=$2
FASTQ_READ2=$3
GENOME=$4

#AFTER TRIMMING (DO NOT EDIT)
TF1=$SAMPLENAME"_1_val_1.fq.gz"
TF2=$SAMPLENAME"_2_val_2.fq.gz"

##EDIT DURING SETUP
###TOOLS:
TRIM=trim_galore
STAR=/mnt/projects/rpd/apps/star-2.5.3a/bin/STAR
HTSEQCOUNT=/mnt/software/bin/htseq-count
SAMTOOLS=/mnt/bin/software

###DATABASE/ANNOTATION:
ADAPTOR=/mnt/projects/wlwtan/cardiac_epigenetics/foolab/jenny/mar2020/rnaseq/analysis_adaptor/illumina.fa
STARIND="/mnt/projects/rpd/genomes/"$GENOME"/star"
STARFASTA="/mnt/projects/rpd/genomes/"$GENOME"/"$GENOME".fa"
GTF="/mnt/projects/rpd/genomes/"$GENOME"/gtf/"$GENOME"_annotation.gtf"

###PARAMETERS:
THREAD=4


##ACTUAL COMMANDS:

mkdir $SAMPLENAME
cd $SAMPLENAME

## Trimming of adaptors and base quality
$TRIM --fastqc --gzip --length 100 --paired $FASTQ_READ1 $FASTQ_READ2


## STAR alignment
$STAR --runThreadN $THREAD --outFileNamePrefix rnaseqtrimmed --outSAMtype BAM Unsorted --genomeDir $STARIND --readFilesCommand zcat --readFilesIn $TF1 $TF2

## sort the bam files by name and count by htseq-count for EdgeR/DESeq analysis
$SAMTOOLS sort -n rnaseqtrimmedAligned.out.bam name_rnaseqtrimmedAligned.out
$HTSEQCOUNT -f bam -r name -s no -m union name_rnaseqtrimmedAligned.out.bam $GTF > count.txt

```  
  
## Requirements / Dependencies:  
- STAR aligner  
- samtools  
- htseq-count  
- trim_galore  
- R (version 3.5 and above)  
- EdgeR  
- DESeq2  
 
  
 
## Part 1: Alignment of RNA-seq reads to Genome
<details>  
<summary> Description of the basic pipeline </summary>  

#### 1. **Trimming of the reads before mapping**  

```
$TRIM --fastqc --gzip --length 100 --paired $FASTQ_READ1 $FASTQ_READ2  

```  

Input files:  
```  
a. FASTQ_READ1 = FASTQ File Read 1  
b. FASTQ_READ2 = FASTQ File Read 2  
```  



  
This command trim the reads:  
a. base quality  
b. adaptors / index  
  
  
   
#### 2. **Mapping of sequencing reads**  
```
$STAR --runThreadN $THREAD --outFileNamePrefix rnaseqtrimmed --outSAMtype BAM Unsorted --genomeDir $STARIND --readFilesCommand zcat --readFilesIn $TF1 $TF2
```  


Input files:  
```
a. TF1 (trimmed fastq file R1)  
b. TF2 (trimmed fastq file R2)  
c. $STARIND = indices required for the mapping of sequencing reads.  

```

The details on various parameters can be found [here](https://github.com/alexdobin/STAR).
  
Running LOG:
```
Mar 13 23:10:09 ..... started STAR run  
Mar 13 23:10:10 ..... loading genome  
Mar 13 23:11:59 ..... started mapping  
Mar 14 02:11:40 ..... finished successfully  
```

  
#### 3. **Generation of RAW gene count using htseq-count**
```
$SAMTOOLS sort -n rnaseqtrimmedAligned.out.bam name_rnaseqtrimmedAligned.out
```   

Input files:
```
a. BAM file (sorted by read name / position )  
b. GTF files  
```

Output file (you need this for most analysis)
```
Filename: count.txt (gene count needed for EdgeR / DESeq)  
Format:  
<Gene ID> <unnormalized count>  
Example:
ENSMUSG00000093369.1	0  
ENSMUSG00000093370.1	3  
ENSMUSG00000093371.1	15  
ENSMUSG00000093373.1	0  
ENSMUSG00000093374.1	0  
__no_feature	9269724  
__ambiguous	97362  
__too_low_aQual	0  
__not_aligned	0  
__alignment_not_unique	3113139  

```
  
Description of the last 5 lines of each files:  
__no_feature: reads (or read pairs) which could not be assigned to any feature   
__ambiguous:  reads (or read pairs) which could have been assigned to more than one feature  
__too_low_aQual: reads (or read pairs) which were skipped due to the -a option  
__not_aligned: reads (or read pairs) in the SAM file without alignment  
__alignment_not_unique: reads (or read pairs) with more than one reported alignment.  
  
  
Running LOG:  
```
...  
31700000 SAM alignment record pairs processed.  
31800000 SAM alignment record pairs processed.  
31900000 SAM alignment record pairs processed.  
32000000 SAM alignment record pairs processed.  
32100000 SAM alignment record pairs processed.  
32200000 SAM alignment record pairs processed.  
32300000 SAM alignment record pairs processed.  
32400000 SAM alignment record pairs processed.  
32500000 SAM alignment record pairs processed.  
32600000 SAM alignment record pairs processed.  
32700000 SAM alignment record pairs processed.  
32800000 SAM alignment record pairs processed.  
32889607 SAM alignment pairs processed.  
```
  
</details>  

<br />  
<br />  
  
## Part 2 Running Differential Expression analysis using edgeR package (using Mouse Sham vs TAC RNA-seq as example)  
  

The differential analysis starts here:  
```
Rscript-3.5.1 edgeR.r <control group> < treatment group>  
```

   
###  Input  
Input file: countfile.txt (hard-coded for now)      
This file is a tab-delimited text file with 3 columns. Header is optional.    
| Column   | Description |
| -------- | ----------- |
| Column 1 | full path to the count file generated by htseq-count |
| Column 2 | sample id                                            |
| Column 3 | condition/grouping information used for DE analysis  |  
  
Example countfile.txt:  
| countfile | sampleid | group |
| --- | --- | --- |
| RMH3869_count.txt | RMH3869 | TAC |
| RMH3871_count.txt | RMH3871 | SHAM |
| RMH3872_count.txt | RMH3872 | TAC |
| RMH3874_count.txt | RMH3874 | SHAM |
| RMH3875_count.txt | RMH3875 | TAC |
| RMH3877_count.txt | RMH3877 | SHAM |

      

### Command  
I usualy run my R script on linux command line. You can also modify the script to make it work on Rstudio or any IDE.  
The command line will require two user arguments for this particular script. These arguments are the sample groupings. In my case, it will be SHAM and TAC.  
I am using R version 3.5.1.  
```
Rscript-3.5.1 single_edgeR.r SHAM TAC
```

<br />
  
### R script (You can copy this into your RStudio if you are not running on linux like me. Dont use single_edgeR.r in that case.)  
This is the skeleton of generic edgeR analysis. You can start with this script first and slowly add/edit according to your experimental designs.  
  

```
filelist = read.table("countfile.txt")
names(filelist)=c("filename","sample","group")
gene2biotype = read.table("/mnt/projects/wlwtan/cardiac_epigenetics/pipeline/rnaseq_pipeline/resource/mm9/mm9_gene2biotype.txt")  
names(gene2biotype)=c("id","hgnc","biotype")  
data = readDGE(file=filelist$filename,group=filelist$group,labels=filelist$sample)  
design = model.matrix(~0+filelist$group)  
filter = apply(data$counts,1,function(x) length(x[x>5])>=2)
filtered = data$counts[ filter, ]
y = DGEList(counts = as.matrix(filtered), group=filelist$group )  
y = calcNormFactors(y)  
y = estimateDisp(y, design)  
fit = glmFit(y,design)
lrt21 = glmLRT(fit,contrast=c(-1,1))  
k21 = as.data.frame(topTags(lrt21,n=10000000))  
k21$hgnc = gene2biotype$hgnc[match(rownames(k21),gene2biotype$id)]  
k21$biotype = gene2biotype$biotype[match(rownames(k21),gene2biotype$id)]  
write.table(k21,file="de_SHAM_TAC.txt",sep="\t",quote=FALSE)  


```
  

### Output of single_edgeR.r   
The R script will generate a few output files.  

1. de_SHAM_TAC.txt  
This is the output of differential gene expression test performed by edgeR. Although there is no definite cut-off to define significance, I use abs(log2FC)>0.5, FDR<0.05 and avg FPKM>1 as a cut-off for significant DEGs. The columns are:  


| Column | Description            |
| :------: | -----------            |
| 1      | GENCODE Gene ID        |
| 2      | log2FC (Group2/Group1) |
| 3 | logCPM across the samples |
| 4 | Test statistics (loglikelihodd) |
| 5 | P-value |
| 6 | FDR |
| 7 | Gene Symbol / HGNC |
| 8 | Biotype |
| 9-whatever | FPKM of all samples |


2. Scatterplot (scatterplot_withtext_SHAM_TAC.png)  
One of the plots that we can look at is scatterplot that represents the correlation between the two conditions (SHAM and TAC). The x-axis in the graph below shows average FPKM of SHAM cardiomyocytes and the y-axis shows average FPKM of TAC cardiomyocytes.  

<p align="center">
  <img height="400" src="https://github.com/lwtan90/wilsonNGS/blob/master/RNAseq/testdata/SHAMTAC/scatterplot_withtext_SHAM_TAC.png">
</p>  




3. Volcano Plot (volcanoplot_withtext_SHAM_TAC.png)  

<p align="center">
  <img height="400" src="https://github.com/lwtan90/wilsonNGS/blob/master/RNAseq/testdata/SHAMTAC/volcanoplot_withtext_SHAM_TAC.png">
</p>  
<br />

4. Unsupervised clustering based on DEGs  
<p align="center">
  <img height="400" src="https://github.com/lwtan90/wilsonNGS/blob/master/RNAseq/testdata/SHAMTAC/degene_heatmap_SHAM_TAC.png">
</p>  
<br />  

5. PCA analysis based on DEGs
<p align="center">
  <img height="400" src="https://github.com/lwtan90/wilsonNGS/blob/master/RNAseq/testdata/SHAMTAC/degenepca_SHAM_TAC.png">
</p>
<br />   

  
  
<br />  
<br />  

## Part 3: Pathway Enrichment Analysis  
### 1. Gene-set Enrichment Analysis (GSEA)  
```
Rscript-3.5.1 fgsea.r de_SHAM_TAC.txt
```  
The analysis was done by using FGSEA R package.  
GSEA analysis takes the logFC / test statistics / pvalue from all the genes in the dataset, and try to find the pathways which are differentially expressed collectively.  

  
Upregulated Pathways:  
  
| pathway | pval | padj | ES | NES |
| :--- | :--: | :--: | :--: | :--: |
| GO_COLLAGEN_TRIMER | 0 | 0.0079 | 0.8466 | 2.682 |
| GO_EXTRACELLULAR_MATRIX_COMPONENT | 0 | 0.0079 | 0.7538 | 2.6497 |
| GO_PROTEINACEOUS_EXTRACELLULAR_MATRIX | 0 | 0.0079 | 0.7161 | 2.6478 |
| GO_EXTRACELLULAR_STRUCTURE_ORGANIZATION | 0 | 0.0079 | 0.7128 | 2.6383 |
| GO_EXTRACELLULAR_MATRIX_BINDING | 0 | 0.0079 | 0.8044 | 2.6162 |
| GO_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENT | 0 | 0.0079 | 0.8255 | 2.6064 |
| GO_EXTRACELLULAR_MATRIX | 0 | 0.0079 | 0.6939 | 2.6011 |
| GO_BASEMENT_MEMBRANE | 0 | 0.0079 | 0.7501 | 2.5799 |
| GO_COLLAGEN_BINDING | 0 | 0.0079 | 0.7807 | 2.564 |
  
  
Downregulated Pathways:  
  

| pathway | pval | padj | ES | NES |
| :--- | :--: | :--: | :--: | :--: |
| GO_MITOCHONDRIAL_ATP_SYNTHESIS_COUPLED_PROTON_TRANSPORT | 0 | 0.0143 | -0.7478 | -2.3421 |
| GO_PROTON_TRANSPORTING_ATP_SYNTHASE_COMPLEX | 0 | 0.0148 | -0.6904 | -2.2669 |
| GO_BRANCHED_CHAIN_AMINO_ACID_METABOLIC_PROCESS | 0 | 0.0148 | -0.6333 | -2.0795 |
| GO_PHOTORECEPTOR_CELL_MAINTENANCE | 0 | 0.016 | -0.5669 | -1.9882 |
| GO_MALE_MEIOSIS | 0 | 0.016 | -0.5459 | -1.9145 |
| GO_NADH_DEHYDROGENASE_ACTIVITY | 0 | 0.0172 | -0.6156 | -2.2153 |
| GO_MITOCHONDRIAL_ELECTRON_TRANSPORT_NADH_TO_UBIQUINONE | 0 | 0.0179 | -0.6089 | -2.2264 |
| GO_NADH_DEHYDROGENASE_COMPLEX | 0 | 0.0179 | -0.5823 | -2.1291 |
| GO_FATTY_ACID_BETA_OXIDATION | 0 | 0.0183 | -0.4836 | -1.8241 |

<br/>
<br/>  
  
### 2. Over-representation Analysis (Enrichment-based) Pathway analysis  
```
Rscript-3.5.1 goseq.r de_SHAM_TAC.txt  
```
The list of differentially expressed genes will be supplied to check for enrichment in a specific geneset or pathway or ontology. The package used here is GOseq. Alternatively, there are also a couple of online free tools that you can use such as David / Enrichr.    

Ontology terms enriched in up-regulated genes:   

| category | over_represented_pvalue | numDEInCat | numInCat | term | ontology |
| :---: | :---: | :---: | :---: | :---: | :---: |
| GO_0030198 | 1.59289631513099e-29 | 58 | 165 | extracellular matrix organization | BP |
| GO_0043062 | 7.25401577545077e-28 | 61 | 198 | extracellular structure organization | BP |
| GO_0007155 | 7.52728027032633e-28 | 130 | 786 | cell adhesion | BP |
| GO_0022610 | 1.89458381420028e-27 | 130 | 794 | biological adhesion | BP |
| GO_0009653 | 1.3895959294514e-21 | 202 | 1808 | anatomical structure morphogenesis | BP |
| GO_0016477 | 4.586619274943e-19 | 127 | 959 | cell migration | BP |
| GO_0001525 | 4.8894400795615e-19 | 72 | 391 | angiogenesis | BP |
| GO_0035295 | 5.83669120066987e-19 | 112 | 791 | tube development | BP |
| GO_0035239 | 9.56339723145197e-19 | 99 | 657 | tube morphogenesis | BP |
| GO_0048870 | 1.80505175084392e-18 | 130 | 1017 | cell motility | BP |


  
 
  
Ontology terms enriched in down-regulated genes:  

| category | over_represented_pvalue | numDEInCat | numInCat | term | ontology |
| :---: | :---: | :---: | :---: | :---: | :---: |
| GO_0044281 | 3.98227858592603e-08 | 71 | 1325 | small molecule metabolic process | BP |
| GO_0006635 | 3.60003697030415e-07 | 11 | 58 | fatty acid beta-oxidation | BP |
| GO_0016054 | 9.25408023716837e-07 | 17 | 154 | organic acid catabolic process | BP |
| GO_0046395 | 9.25408023716837e-07 | 17 | 154 | carboxylic acid catabolic process | BP |
| GO_0009062 | 1.15898676473624e-06 | 12 | 78 | fatty acid catabolic process | BP |
| GO_0008016 | 1.25401426261223e-06 | 16 | 137 | regulation of heart contraction | BP |
| GO_0043436 | 2.97533226290676e-06 | 40 | 668 | oxoacid metabolic process | BP |
| GO_0006082 | 3.93851000867858e-06 | 40 | 676 | organic acid metabolic process | BP |
| GO_0044242 | 4.12629854308222e-06 | 15 | 136 | cellular lipid catabolic process | BP |
| GO_0072329 | 4.19219279151121e-06 | 12 | 88 | monocarboxylic acid catabolic process | BP |


  

