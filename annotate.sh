#!/bin/sh
#https://github.com/TransDecoder/TransDecoder/wiki
module load Miniconda3
source activate transdecoder
#TransDecoder.LongOrfs -t Trinity.fasta
TransDecoder.Predict -t Trinity.fasta



makeblastdb -in Trinity.fasta.transdecoder.pep -dbtype prot
blastx -query  12genelucila.fasta  -db Trinity.fasta.transdecoder.pep -outfmt 7   >  blast_results