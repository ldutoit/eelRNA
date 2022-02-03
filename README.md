# README

This project is associated to:

Lucila Babio, L., Lokman, P. M., Damsteegt, E. L., & Dutoit, L. Are Cell Junctions Implicated in the Regulation of Vitellogenin
Uptake? Insights from an RNAseq-based Study in Eel, Anguilla australis. Cells, In Press


## Summary

We build the transcriptome of the NZ eel denovo using Trinity, we then performed a quick annotation using Trinotate. 

We then run differential expression between 6 Yellow and 6 silver eels with 100bp paired-end data. A

ll the differentially expressed genes with statistics and counts are saved in [results_files/DE_results.txt](results_files/DE_results.txt). Statistics and counts for ALL Genes are in [results_files/ALLgenes_results.txt](results_files/ALLgenes_results.txt).

The Gene ontology enrichment analysis is saved in [results_files/diff.GOseq.enriched](results_files/diff.GOseq.enriched) and [diff.GOseq.depleted](results_files/diff.GOseq.depleted). 

Finally, all the transcripts fasta sequences are in [results_files/Trinity.fasta.gz](results_files/Trinity.fasta.gz)

## Description

The raw data is saved on the SRA at project PRJNA785278.

## Analyses

### Building the transcriptome and obtaining counts.

The transcriptome is built in [trinityrun.md](trinityrun.md) all the way to obtaining gene-level counts for the 12 samples. Annotation of the transcriptome was done using [annotate.md](annotate.md). That creates the Trinotate output [denovoassembly_annotation_report.xls](denovoassembly_annotation_report.xls) Which is the complete annotation at the transcript level. 


### Differential Expression (DE) Analysis

The differential expression analysis is done in [DE_analysis_eels.md](DE_analysis_eels.md). 

### GO Enrichment

We  extracted GO terms at the gene level  and performed the GO enrichment analysis in the file [GOenrichment.md](GOenrichment.md). The outputs are in the text files [results_files/diff.GO.enriched](diff.GO.enriched) and  [results_files/diff.GO.depleted](diff.GO.depleted).

