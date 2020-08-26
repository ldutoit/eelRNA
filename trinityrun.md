# Trinity run

Trinity is run on the mahuika server on the nobackup partition nobackup following nesi instructions https://support.nesi.org.nz/hc/en-gb/articles/360000980375-Trinity


```
cd /nesi/nobackup/uoo00116/eelrna
```


After checking the quality of the data within

using Fastqc.

I learn that the data is very good with 35Mio paired end reads per sample.

ln -s ~/projects/eelRNA/source_files/ .
cd /nesi/nobackup/uoo00116/eelrna/source_files/CleandatafromBGI

I create a sample files  as per below (samples_file.txt)


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



Then I run trinity with the two strand option given I am not sure which one I should use from the BGI report.


```
#!/bin/bash -e
#SBATCH --job-name=trinity-phase1_RF
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
  --samples_file samples_file.txt --SS_lib_type RF   --output RF_trinity_output
```



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



### Obtaining iscount matrix


The next step is to obtain a count matrix
I used the tutorial at ...
https://southgreenplatform.github.io/trainings/trinityTrinotate/TP-trinity/



# create a salmon_outdir 

TO Adapt

# salmon
util/align_and_estimate_abundance.pl \
--transcripts FR_trinity_output/Trinity.fasta \
--seqType fq \
--samples_file samples_file.txt \
--est_method salmon \
--trinity_mode \
--prep_reference > salmon_align_and_estimate_abundance.log 2>&1 &


$path_to_trinity/util/abundance_estimates_to_matrix.pl \
--est_method salmon \
--out_prefix Trinity_trans \
--name_sample_by_basedir \
--gene_trans_map none \
CENPK_rep1/quant.sf \
CENPK_rep2/quant.sf \
CENPK_rep3/quant.sf \
Batch_rep1/quant.sf \
Batch_rep2/quant.sf \
Batch_rep3/quant.sf 


Buscp