
First trun transdecoder
```
#!/bin/sh
#https://github.com/TransDecoder/TransDecoder/wiki
module load Miniconda3
source activate transdecoder
#TransDecoder.LongOrfs -t Trinity.fasta
TransDecoder.Predict -t Trinity.fasta
```


Annotate using trinotate

#https://southgreenplatform.github.io/trainings//files/AA-SG-ABiMS2019/RNASeq_denovo_Montpellier_092019_3_annotation.pdf

```
module load Miniconda3source activate transdecoder # self install

Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate
Trinotate Trinotate.sqlite init --gene_trans_map Trinity.fasta.gene_trans_map \
--transcript_fasta Trinity.fasta --transdecoder_pep Trinity.fasta.transdecoder.pep
```


Do the blasting:

```
#Blast
module load Blast
makeblastdb --in  uniprot_sprot.pep --db-type prot
blastx -query Trinity.fasta -db uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastx.outfmt6
blastp -query Trinity.fasta.transdecoder.pep -db uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastp.outfmt6

# hmm
gunzip -d Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
hmmscan --cpu 4 --domtblout TrinotatePFAM.out Pfam-A.hmm \
Trinity.fasta.transdecoder.pep > pfam.log #Pfam come from Build_Trinotate_Boilerplate_SQLite_db above


#signalp weird installation I need to fix the paths and use two symlinks

ln -s /home/ludovic.dutoit/repos/softwares/signalp/signalp-5.0b/bin/ .

bin/signalp -format short -prefix signalp -fasta  Trinity.fasta.transdecoder.pep
```
	Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp.outfmt6 
```