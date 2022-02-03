# Trinity run

After checking the quality of the data within using Fastqc, I confirm that the data is very good with 35Mio paired end reads per sample.

```
ln -s ~/projects/eelRNA/source_files/ .
cd /nesi/nobackup/uoo00116/eelrna/source_files/CleandatafromBGI
```

I create a sample file as per below, saved in  [samples_file.txt](samples_file.txt)


```
silver	silver_rep1	source_files/CleandatafromBGI/1LS_1.fq.gz	source_files/CleandatafromBGI/1LS_2.fq.gz
yellow	yellow_rep1	source_files/CleandatafromBGI/1LY_1.fq.gz	source_files/CleandatafromBGI/1LY_2.fq.gz
silver	silver_rep2	source_files/CleandatafromBGI/2LS_1.fq.gz	source_files/CleandatafromBGI/2LS_2.fq.gz
yellow	yellow_rep2	source_files/CleandatafromBGI/2LY_1.fq.gz	source_files/CleandatafromBGI/2LY_2.fq.gz
silver	silver_rep3	source_files/CleandatafromBGI/3LS_1.fq.gz	source_files/CleandatafromBGI/3LS_2.fq.gz
yellow	yellow_rep3	source_files/CleandatafromBGI/3LY_1.fq.gz	source_files/CleandatafromBGI/3LY_2.fq.gz
silver	silver_rep4	source_files/CleandatafromBGI/4LS_1.fq.gz	source_files/CleandatafromBGI/4LS_2.fq.gz
yellow	yellow_rep4	source_files/CleandatafromBGI/4LY_1.fq.gz	source_files/CleandatafromBGI/4LY_2.fq.gz
silver	silver_rep5	source_files/CleandatafromBGI/5LS_1.fq.gz	source_files/CleandatafromBGI/5LS_2.fq.gz
yellow	yellow_rep5	source_files/CleandatafromBGI/5LY_1.fq.gz	source_files/CleandatafromBGI/5LY_2.fq.gz
silver	silver_rep6	source_files/CleandatafromBGI/6LS_1.fq.gz	source_files/CleandatafromBGI/6LS_2.fq.gz
yellow	yellow_rep6	source_files/CleandatafromBGI/6LY_1.fq.gz	source_files/CleandatafromBGI/6LY_2.fq.gz
```

Then I run trinity in 2 phases as per my specific [cluster instructions](https://support.nesi.org.nz/hc/en-gb/articles/360000980375-Trinity). This is simply a locally optimised way of running the main `Trinit` denovo assembly command.


```
#!/bin/bash -e
#SBATCH --job-name=trinity-phase1_FR
#SBATCH --account=uoo00116   # your NeSI project code
#SBATCH --time=30:00:00       # maximum run time
#SBATCH --ntasks=1            # always 1
#SBATCH --cpus-per-task=16    # number of threads to use for Trinity
#SBATCH --mem=220G            # maximum memory available to Trinity
#SBATCH --partition=bigmem    # based on memory requirements
#SBATCH --hint=nomultithread  # disable hyper-threading

# load a Trinity module
module load Trinity/2.8.5-gimkl-2018b

# run trinity, stop before phase 2
srun Trinity --no_distributed_trinity_exec \
  --CPU ${SLURM_CPUS_PER_TASK} --max_memory 200G \
  --samples_file samples_file.txt --SS_lib_type FR   --output FR_trinity_output
```

And phase2:
```
#!/bin/bash -e
#SBATCH --job-name=trinity-phase2grid
#SBATCH --account=uoo00116  # your NeSI project code
#SBATCH --time=30:00:00      # enough time for all sub-jobs to complete
#SBATCH --ntasks=1           # always 1 - this is the master process
#SBATCH --cpus-per-task=1    # always 1
#SBATCH --mem=20G            # memory requirements for master process
#SBATCH --partition=hugemem   # submit to an appropriate partition
#SBATCH --hint=nomultithread

# load Trinity and HPC GridRunner
module load Trinity/2.8.5-gimkl-2018b
module load HpcGridRunner/20181005

# run Trinity - this will be the master HPC GridRunner process that handles
#   submitting sub-jobs (batches of commands) to the Slurm queue
srun Trinity --CPU ${SLURM_CPUS_PER_TASK} --max_memory 20G \
  --grid_exec "hpc_cmds_GridRunner.pl --grid_conf ${SLURM_SUBMIT_DIR}/SLURM.conf -c" \
   --seqType fq --samples_file samples_file.txt --SS_lib_type FR   --output FR_trinity_output
```

### Obtaining iscount matrix


The next step is to obtain a count matrix
I used the tutorial at:
https://southgreenplatform.github.io/trainings/trinityTrinotate/TP-trinity/

The combination of modules is a bit over the top (due to compiler issues) but it works for me.

```
#user specific environment
module purge
module load Trinity/2.8.4-gimkl-2017a
module load Miniconda3
source activate transdecoder #modu that contains rm #rsem and other of my trinity utils
module load  Trinity/2.8.4-gimkl-2017a SAMtools/1.8-gimkl-2017a Bowtie/1.2.0-gimkl-2017a RSEM/1.3.1-gimkl-2017a

#Trinity command
align_and_estimate_abundance.pl \
--transcripts FR_trinity_output/Trinity.fasta \
--seqType fq \
--samples_file samples_file.txt \
--est_method RSEM --aln_method bowtie2 \
--trinity_mode \
--prep_reference \
--thread_count 16 \
--coordsort_bam > bowtie-rsem_align_and_estimate_abundance.log 
```

We obtain gene counts organised in 12 folder by samples names. We group them in one folder ysing links:

```bash
mkdir RSEM_results
cd RSEM_results
ln -s  ../silver_rep1/RSEM.genes.results silver_rep1.txt
ln -s  ../silver_rep2/RSEM.genes.results silver_rep2.txt
ln -s  ../silver_rep3/RSEM.genes.results silver_rep3.txt
ln -s  ../silver_rep4/RSEM.genes.results silver_rep4.txt
ln -s  ../silver_rep5/RSEM.genes.results silver_rep5.txt
ln -s  ../silver_rep6/RSEM.genes.results silver_rep6.txt
ln -s  ../yellow_rep1/RSEM.genes.results yellow_rep1.txt
ln -s  ../yellow_rep2/RSEM.genes.results yellow_rep2.txt
ln -s  ../yellow_rep3/RSEM.genes.results yellow_rep3.txt
ln -s  ../yellow_rep4/RSEM.genes.results yellow_rep4.txt
ln -s  ../yellow_rep5/RSEM.genes.results yellow_rep5.txt
ln -s ../yellow_rep6/RSEM.genes.results yellow_rep6.txt
```			

In order to create a singleclean matrix we go into R and use the tximport module

from within the directory  ```gene_counts/RSEM_results```

```R
library("tximport")
files <- dir()
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
head(txi.rsem$counts)
colnames(txi.rsem$counts) <- gsub(".txt","",files)
head(txi.rsem$counts)
#I checked manually that the colnames match the counts in the single files

write.table(txi.rsem$counts,"RSEM_gene_counts.txt",row.names=T,col.names=T,sep="\t")
```


I saved the count matrix in this repository [results_files/RSEM_gene_counts.txt](results_files/RSEM_gene_counts.txt)

