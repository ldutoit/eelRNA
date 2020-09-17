
First trun transdecoder

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

# hmm
gunzip -d Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
hmmscan --cpu 4 --domtblout TrinotatePFAM.out Pfam-A.hmm \
Trinity.fasta.transdecoder.pep > pfgit#Pfam come from Build_Trinotate_Boilerplate_SQLite_db above
```


```bash
Trinotate Trinotate.sqlite init --gene_trans_map Trinity.fasta.gene_trans_map \
--transcript_fasta Trinity.fasta --transdecoder_pep trinity.fasta.transdecoder.pep

Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp.outfmt6
Trinotate Trinotate.sqlite LOAD_swissprot_blastx
Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite report > denovoassembly_annotation_report.xls
```

That is the final assembly report.

I want to add two columns - one that says wether it is expressed passing filter in the differential expression analsyis - one whether it is DE


I'll start by adding these two columns:

in bash, I extract the second column of the annotation which is all the isoforms.I do that so I don't have to worry about the formatting of all the other columns which won't fit into R.

```bash
 cut -f 2 denovoassembly_annotation_report.xls > tempisoformsinorder.txt
```

Then I use R to generate a file that can be pasted to the existing annotation.


```
load("countsand_logFCforplotting.RData") # load results 

gene_names<-read.table("tempisoformsinorder.txt",h=T)
passing_expressionDE_filter = rep("NO",times=length(gene_names[,1]))
padj =rep("NO",times=length(gene_names[,1]))

for (transcript_nr in 1:length(gene_names[,1])){
  if (transcript_nr%%100==0){print(transcript_nr)}
  if (as.character(gene_names[transcript_nr,1])%in%as.character(rownames(results))){
    passing_expressionDE_filter[transcript_nr]<-"YES"
    padj[transcript_nr]<-results$adj.P.Val[which(as.character(gene_names[transcript_nr,1])==as.character(rownames(results)))]
  }
  }

write.table(cbind(passing_expressionDE_filter,padj),"tempoutput",sep="\t",quote=F)


```

I paste both files:

```
paste denovoassembly_annotation_report.xls tempoutput >  denovoassembly_annotation_report_withDE.xls
#rm  tempisoformsinorder tempoutput
```

