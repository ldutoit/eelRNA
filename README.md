# README

This project is associated to:

Lucila Babio, L., Lokman, P. M., Damsteegt, E. L., & Dutoit, L. Are Cell Junctions Implicated in the Regulation of Vitellogenin
Uptake? Insights from an RNAseq-based Study in Eel, Anguilla australis. Cells, In Press

Feel free to ask for help to dutoit.ludovic@gmail.com when accessing this repository. 

## Summary

We build the transcriptome of the NZ eel denovo using Trinity, we then performed a quick annotation using Trinotate. 

We then run differential expression between 6 Yellow and 6 silver eels with 100bp paired-end data. 

All the differentially expressed genes with statistics and counts are saved in [results_files/DE_results.txt](results_files/DE_results.txt). Statistics and counts for ALL Genes are in [results_files/ALLgenes_results.txt](results_files/ALLgenes_results.txt).

The Gene ontology enrichment analysis is saved in [results_files/diff.GOseq.enriched](results_files/diff.GOseq.enriched) and [diff.GOseq.depleted](results_files/diff.GOseq.depleted). 

Finally, all the transcripts fasta sequences are in [results_files/Trinity.fasta.gz](results_files/Trinity.fasta.gz)

## Description

The raw data is saved on the SRA at project PRJNA785278. The unattoated transcriptome is [results_files/Trinity.fasta.gz](results_files/Trinity.fasta.gz). The annotation Trinotate output is [results_files/denovoassembly_annotation_report.xls.gz](results_files/denovoassembly_annotation_report.xls.gz). 


## Analyses

For detailed methods, refer to the published article above.

### Building the transcriptome and obtaining counts.

The transcriptome is built in [trinityrun.md](trinityrun.md) all the way to obtaining gene-level counts for the 12 samples. Annotation of the transcriptome was done using [annotate.md](annotate.md). That creates the Trinotate output [results_files/denovoassembly_annotation_report.xls](results_files/denovoassembly_annotation_report.xls.gz) Which is the complete annotation at the transcript level. 


### Differential Expression (DE) Analysis

All necessary files are in this repository. The differential expression analysis is done in [DE_analysis_eels.md](DE_analysis_eels.md). 

### GO Enrichment

We  extracted GO terms at the gene level  and performed the GO enrichment analysis in the file [GOenrichment.md](GOenrichment.md). The outputs are in the text files [results_files/diff.GO.enriched](allDE.GOseq.enriched) and  [results_files/diff.GO.depleted](allDE.GOseq.depleted).

