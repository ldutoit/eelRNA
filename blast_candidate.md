This short commands shows how to blast candidate genes of interest to the transcriptome using predicted peptide sequence obtained using transdecoder in [annotate.md](annotate.md)

```
makeblastdb -in Trinity.fasta.transdecoder.pep -dbtype prot
blastx -query  12genelucila.fasta  -db Trinity.fasta.transdecoder.pep -outfmt 7   >  blast_results.txt
```

makeblastdb -in Trinity.fasta.transdecoder.pep -dbtype prot
blastx -query  12genelucila.fasta  -db Trinity.fasta.transdecoder.pep -outfmt 7  -max_target_seqs 2  >  blast_results.txt

