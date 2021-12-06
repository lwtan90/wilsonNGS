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


