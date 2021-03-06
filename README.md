# README

This project aim at comparing the transcriptome of yellow and silver NZ eels.

## Summary

We build the transcriptome of the NZ eel denovo using Trinity. We then run differential expression between 6 Yellow and 6 silver eels with 100bp paired-end data. All the differentially expressed genes with statistics and counts are saved in [results_files/DE_results.txt](results_files/DE_results.txt). Statistics and counts for ALL Genes are in [results_files/ALLgenes_results.txt](results_files/ALLgenes_results.txt).

The Gene ontology enrichment analysis is saved in [results_files/diff.GOseq.enriched](results_files/diff.GOseq.enriched) and [diff.GOseq.depleted](results_files/diff.GOseq.depleted). They are output of GOSeq and one should pay attention to the False Discovery Rate values!

Finally, all the transcripts fasta sequences are in [results_files/Trinity.fasta.gz](results_files/Trinity.fasta.gz)

## Description

The raw data is saved on the HCS of otago university: 
```
/sci-bioinformatics-project-archive/Transcriptome data Lucila 2019
```
contact dutoit.ludovic@gmail.com for access

## Analyses

### Building the transcriptome and obtaining counts

The transcriptome is built in [trinityrun.md](trinityrun.md) all the way to obtaining gene-level counts for the 12 samples.

### Differential Expression (DE) Analysis

The differential expression analysis is done in [DE_analysis_eels.md](DE_analysis_eels.md). 

### GO Enrichment

Annotation of the transcriptome is done in [annotate.md](annotate.md). That creates the file denovoassembly_annotation_report.xls Which is the complete annotation at the transcript level. 

We then extract GO terms at the gene level  and performed the GO enrichment analysis in the file [GOenrichment.md](GOenrichment.md)



## Utilities

A couple of utilities to be able to dig in the data using R.



[plotting_genes.R](plotting_genes.R) allows to plot the transcripts read counts. It needs both ggplot2 (the first time you may need to install it using install.packages(“tidyverse”) and the file countsand_logFCforplotting.RData inside the same folder than the script (that file contains the results of the analysis and the counts table as R objects).

[plotting_transcript.R](plotting_transcript.R)The same as plotting genes but utilising salmon transcript-level counts.


[extract_sequences.R](extract_sequences.R) This script extract ssequences from the transcript files or the peptide files and save it one sequence file. It works with the file Trinity.fasta ( mRNAs) and the file Trinity.fasta.transdecoder.pep ( contact author for access).


[blast_candidate.md](blast_candidate.md) allows to blast genes using an example of 12 genes from Lucila.
