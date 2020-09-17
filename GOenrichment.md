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

extract gene lengths to be able to account for gene length in GOseq annotation.


```
/opt/nesi/CS400_centos7_bdw/Trinity/2.8.5-gimkl-2018b/trinityrnaseq-Trinity-v2.8.5/util/misc/fasta_seq_length.pl Trinity.fasta > Trinity.fasta.seq_lens
```
```
/opt/nesi/CS400_centos7_bdw/Trinity/2.8.5-gimkl-2018b/trinityrnaseq-Trinity-v2.8.5/util/misc/TPM_weighted_gene_length.py  \
         --gene_trans_map trinity_out_dir/Trinity.fasta.gene_trans_map \
         --trans_lengths Trinity.fasta.seq_lens \
         --TPM_matrix isoforms.TMM.EXPR.matrix > Trinity.gene_lengths.txt
```
