# Annotate

### First trun transdecoder

```bash
#!/bin/sh
#https://github.com/TransDecoder/TransDecoder/wiki
module load Miniconda3
source activate transdecoder
#TransDecoder.LongOrfs -t Trinity.fasta
TransDecoder.Predict -t Trinity.fasta
```


### Annotate using trinotate

https://southgreenplatform.github.io/trainings//files/AA-SG-ABiMS2019/RNASeq_denovo_Montpellier_092019_3_annotation.pdf

```bash
module load Miniconda3source activate transdecoder # self install

Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate
Trinotate Trinotate.sqlite init --gene_trans_map Trinity.fasta.gene_trans_map \
--transcript_fasta Trinity.fasta --transdecoder_pep Trinity.fasta.transdecoder.pep
```


Do the blasting:

```bash
#Blast
module load Blast
makeblastdb --in  uniprot_sprot.pep --db-type prot
blastx -query Trinity.fasta -db uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastx.outfmt6
blastp -query Trinity.fasta.transdecoder.pep -db uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastp.outfmt6


```


```bash
Trinotate Trinotate.sqlite init --gene_trans_map Trinity.fasta.gene_trans_map \
--transcript_fasta Trinity.fasta --transdecoder_pep trinity.fasta.transdecoder.pep

Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp.outfmt6
Trinotate Trinotate.sqlite LOAD_swissprot_blastx
Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite report > denovoassembly_annotation_report.xls
```

That is the final assembly report.


