# README

This project aim at comparing the transcriptome of yellow and silver NZ eels.

## Summary

We build the transcriptome of the NZ eel denovo using Trinity. We then run differential expression between 6 Yellow and 6 silver eels with 100bp paired-end data.

## Description
 of Bruce Robertson. The SNP were called by Yasmin Foster.

The raw data is saved on the HCS of otago university: 

```
/sci-bioinformatics-project-archive/Transcriptome data Lucila 2019
```
contact dutoit.ludovic@gmail.com for access

## Analyses

The transcriptome is built in trinityrun.md



## Utilities

A couple of utilities to be able to dig in the data using R.

[plotting_transcript.R](plotting_transcript.R) allows to plot the transcripts read counts. It needs both ggplot2 (the first time you may need to install it using install.packages(“tidyverse”) and the file countsand_logFCforplotting.RData inside the same folder than the script (that file contains the results of the analysis and the counts table as R objects).

[extract_sequences.R](extract_sequences.R) This script extract ssequences from the transcript files or the peptide files and save it one sequence file. It works with the file Trinity.fasta ( mRNAs) and the file Trinity.fasta.transdecoder.pep ( contact author for access).
