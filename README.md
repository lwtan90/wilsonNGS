
# Wilson's Bioinformatics Notebook  
Started in March 2020

## RNA-seq  
### Description  
The RNA-seq pipeline takes in paired-end FASTQ file (gzipped / unzipped), and produced gene count, differential gene expression analysis results, isoform estimation, splicing analysis, promoter motif analysis and pathway analysis.
The pipeline was designed to run on SGE cluster.  

### Dependencies  
1. trim_galore  
2. STAR  
3. Kallisto  
4. RSEM  
5. htseq-count  
6. R  
7. EdgeR  
8. DESeq  
9. Sleuth  
