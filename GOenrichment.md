### go enchrichment analysis


[results_files/denovoassembly_annotation_report.xls.gz]([results_files/denovoassembly_annotation_report.xls.gz) is the original trinotate output. Note that of the code below is specific to my installation of trinotate.

https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-GOSeq contains the methodology fr the go enrichment analysis and was used as a template.



### Extract go terms

```
module load Miniconda3 #specific to my installation of trinotate
source activate transdecoder  #specific to my installation of trinotate
extract_GO_assignments_from_Trinotate_xls.pl --gene \
                         --Trinotate_xls  denovoassembly_annotation_report.xls \
                         -G --include_ancestral_terms \
                         > go_annotations.txt
````

##Extract GENE lengths files

extract gene lengths to be able to account for gene length in GOseq annotation.

```
/opt/nesi/CS400_centos7_bdw/Trinity/2.8.5-gimkl-2018b/trinityrnaseq-Trinity-v2.8.5/util/misc/fasta_seq_length.pl FR_trinity_output/Trinity.fasta > Trinity.fasta.seq_lens
```
```

```

 isoforms.TMM.EXPR.matrix  is created using counts generated using salmon and the script:


```

/opt/nesi/CS400_centos7_bdw/Trinity/2.8.5-gimkl-2018b/trinityrnaseq-Trinity-v2.8.5/util/align_and_estimate_abundance.pl \
--transcripts FR_trinity_output/Trinity.fasta \
--seqType fq \
--samples_file samples_file.txt \
--est_method salmon \
--trinity_mode \
--prep_reference \
 > salmon_align_and_estimate_abundance.log 2>&1 &


/opt/nesi/CS400_centos7_bdw/Trinity/2.8.5-gimkl-2018b/trinityrnaseq-Trinity-v2.8.5/util/abundance_estimates_to_matrix.pl \
--est_method salmon \
--out_prefix Trinity_trans \
--name_sample_by_basedir \
--gene_trans_map none \
silver_rep1/quant.sf \
silver_rep2/quant.sf \
silver_rep3/quant.sf \
silver_rep4/quant.sf \
silver_rep5/quant.sf \
silver_rep6/quant.sf \
yellow_rep1/quant.sf \
yellow_rep2/quant.sf \
yellow_rep3/quant.sf \
yellow_rep4/quant.sf \
yellow_rep5/quant.sf \
yellow_rep6/quant.sf 
```

Finish the gene lengths:

```
/opt/nesi/CS400_centos7_bdw/Trinity/2.8.5-gimkl-2018b/trinityrnaseq-Trinity-v2.8.5/util/misc/TPM_weighted_gene_length.py  \
         --gene_trans_map FR_trinity_output/Trinity.fasta.gene_trans_map \
         --trans_lengths FR_trinity_output/Trinity.fasta.seq_lens \
         --TPM_matrix Trinity_trans.isoform.TMM.EXPR.matrix >Trinity.gene_lengths.txt
```

#I create the factor labelling out of the DE analysis genes in the github repos:

```
 cut -f 2  ~/repos/scripts/eelRNA/results_files/DE_results.txt  | tail -n 10348  |  awk '{print "diff\t",$0}' | sed -s 's/\"//g'    >factor_labeling.txt
```


# #Install GOSeq dependencies ( not covered in tutorial)


 I installed goseq using conda inside the transdecoder conda environment:

```bash
 conda install -c bioconda bioconductor-goseq 

 module load R/3.4.2-gimkl-2017a
```

Then I installed qvalue inside R using:

```r
source("http://bioconductor.org/biocLite.R")
biocLite("qvalue") # that on
biocLite("goseq")  
```

## RUN THE GO ANALYSIS:

```bash
/opt/nesi/CS400_centos7_bdw/Trinity/2.8.5-gimkl-2018b/trinityrnaseq-Trinity-v2.8.5/Analysis/DifferentialExpression/run_GOseq.pl \
                       --factor_labeling  factor_labeling.txt \
                       --GO_assignments go_annotations.txt \
                       --lengths Trinity.gene_lengths.txt \
                       --background  backgroundGO.txt
mv diff.GOseq.depleted  ~/repos/scripts/eelRNA/results_files/allDE.GOseq.depleted 
mv diff.GOseq.enriched  ~/repos/scripts/eelRNA/results_files/allDE.GOseq.enriched                 
```



**Upregulated only**

```bash
cat ~/repos/scripts/eelRNA/upregulated.txt |  awk '{print "diff\t",$0}'    >factor_labelingupregulated.txt
/opt/nesi/CS400_centos7_bdw/Trinity/2.8.5-gimkl-2018b/trinityrnaseq-Trinity-v2.8.5/Analysis/DifferentialExpression/run_GOseq.pl \
                       --factor_labeling  factor_labelingupregulated.txt \
                       --GO_assignments go_annotations.txt \
                       --lengths Trinity.gene_lengths.txt \
                       --background  backgroundGO.txt
mv diff.GOseq.depleted  ~/repos/scripts/eelRNA/results_files/logFCbiggerthan0_diff.GOseq.depleted 
mv diff.GOseq.enriched  ~/repos/scripts/eelRNA/results_files/logFCbiggerthan0_diff.GOseq.enriched

```
**down regulated only**

```bash
cat ~/repos/scripts/eelRNA/downregulated.txt |  awk '{print "diff\t",$0}'    >factor_labelingdownregulated.txt
/opt/nesi/CS400_centos7_bdw/Trinity/2.8.5-gimkl-2018b/trinityrnaseq-Trinity-v2.8.5/Analysis/DifferentialExpression/run_GOseq.pl \
                       --factor_labeling  factor_labelingdownregulated.txt \
                       --GO_assignments go_annotations.txt \
                       --lengths Trinity.gene_lengths.txt \
                       --background  backgroundGO.txt                       
mv diff.GOseq.depleted  ~/repos/scripts/eelRNA/results_files/logFClowerthan0_diff.GOseq.depleted 
mv diff.GOseq.enriched  ~/repos/scripts/eelRNA/results_files/logFClowerthan0_diff.GOseq.enriched
```
I save those files in this repository as [results_files/allDE.GOseq.depleted](results_files/allDE.GOseq.depleted) and [results_files/dallDE.GOseq.depleted](results_files/allDE.GOseq.enriched)  
