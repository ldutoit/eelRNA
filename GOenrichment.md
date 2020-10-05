### go enchrichment analysis


denovoassembly_annotation_report.xls is the original trinotate output

https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-GOSeq contains the methodology fr the go enrichment analysis



### Extract go terms

```
module load Miniconda3
source activate transdecoder
extract_GO_assignments_from_Trinotate_xls.pl --gene \
                         --Trinotate_xls  denovoassembly_annotation_report.xls \
                         -G --include_ancestral_terms \
                         > go_annotations.txt
````

##Extract GENE lengths files

extract gene lengths to be able to account for gene length in GOseq annotation.

```
/opt/nesi/CS400_centos7_bdw/Trinity/2.8.5-gimkl-2018b/trinityrnaseq-Trinity-v2.8.5/util/misc/fasta_seq_length.pl FR_trinity_output/Trinity.fasta > Trinity.fasta.seq_lens
```
```
/opt/nesi/CS400_centos7_bdw/Trinity/2.8.5-gimkl-2018b/trinityrnaseq-Trinity-v2.8.5/util/misc/TPM_weighted_gene_length.py  \
         --gene_trans_map FR_trinity_output/Trinity.fasta.gene_trans_map \
         --trans_lengths FR_trinity_output/Trinity.fasta.seq_lens \
         --TPM_matrix Trinity_trans.isoform.TPM.not_cross_norm > !! THIS TPM MATRIX MIGHT BE WRONG and Trinity.fasta.gene_trans_map too Trinity.gene_lengths.txt
```

I create the factor labelling out of the DE analysis genes in the github repos:


```
 cut -f 1  ~/repos/scripts/eelRNA/DE_results.txt  | tail -n 10348  |  awk '{print "diff\t",$0}' -   >factor_labeling.txt
```

 NOTE 10348 because that is the number of DE genes, there is 10349 lines with a header 

##Install GOSeq dependencies ( not covered in tutorial)


 I installed goseq using conda inside the transdecoder conda environment:

```bash
 conda install -c bioconda bioconductor-goseq 
```

Then I installed qvalue inside R using:

```r
source("http://bioconductor.org/biocLite.R")
biocLite("qvalue") # that one 
```

## RUN THE GO ANALYSIS:

```bash
/opt/nesi/CS400_centos7_bdw/Trinity/2.8.5-gimkl-2018b/trinityrnaseq-Trinity-v2.8.5/Analysis/DifferentialExpression/run_GOseq.pl \
                       --factor_labeling  factor_labeling.txt \
                       --GO_assignments go_annotations.txt \
                       --lengths Trinity.gene_lengths.txt \
                       --background  backgroundGO.txt
```
!!!I NEED TO CLEAN THOSE GENE LENGTHS BASED ON LINE 31 of this file@@@
!!!I NEED TO CHECK WHAT IS THE UNIVERSE/BACKGROUND FOR DE GENES
I save those files in this repository as [diff.GOseq.depleted](diff.GOseq.depleted) and [diff.GOseq.enriched](diff.GOseq.enriched)

