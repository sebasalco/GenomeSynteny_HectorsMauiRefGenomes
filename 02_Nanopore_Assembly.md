# Nanopore long-reads assembly pipeline. Each step was performed independently for each of the Hector's and Māui dolphins assemblies. The scripts won't have the complete path to files and will be focused on Hector's files to exemplify their use.
---
## 1. Guppy basecalling
We used Guppy version 6.4.6 to basecall all the fast5 files produced by minions. It was run in GPU partition. 

`Script for Guppy`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  4
#SBATCH --job-name       GuppyBaseCalling
#SBATCH --mem            20G
#SBATCH --time           10:00:00
#SBATCH --gpus-per-node  A100:1
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load ont-guppy-gpu/6.4.6

# 1. Run Guppy
guppy_basecaller -i /path_to_fast5_files_directory/fast5 \ #Path to the folder containing all the fast5 files from Nanopore
		-s /path_for_output/GuppyHectors \ #Path to the desired output directory
		--config /path_to_configuration_file/dna_r9.4.1_450bps_hac.cfg  \ #Flowcell and kit configuration FLO-MIN106, SQK-LSK110
		--device auto \ #Prioritize GPU
		--recursive \ #Look for all fast5 files recursively
		--num_callers 4 -x auto \ #Number of parallel basecallers
		--trim_barcodes \ #Trim the barcodes from the output sequence in the FastQ files
		--disable_qscore_filtering  #Disable qscore filtering in pass and fail directories

# 2. Concatenate all fastq files into a single file

cat /path_for_output/GuppyHectors/fastq/*.fastq > all_guppyHectors.fastq
```

## 2. PycoQC
Guppy output generates a summary file 'sequencing_summary.txt'. We use that file to analyse the quality of the base calling before proceeding to the next steps

`Script for PycoQC`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  1
#SBATCH --job-name       PycoQC
#SBATCH --mem            10G
#SBATCH --time           01:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load pycoQC/2.5.2-gimkl-2020a-Python-3.8.2

pycoQC	--summary_file /path_to_guppy_sequencing_summary_output/sequencing_summary.txt 
	--html_outfile /path_for_output/pycoQC_output.html
```
## 3. Nanolyse
Remove lambda DNA from our fastq files.

`Script for Nanolyse`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  1
#SBATCH --job-name       Nanolyse
#SBATCH --mem            15G
#SBATCH --time           05:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load NanoLyse/1.2.0-gimkl-2020a

cat /path_to_allfastq_file/all_guppyHectors.fastq | NanoLyse --reference /path_to_reference_lambda_genome/lambda_3.6kb.fasta | gzip > hectors_nanopore_merged_filtered.fastq.gz #Remove the lambda DNA using the lambda reference genome and gzip the results in a new file
```
## 4. Porechop
We remove any adapter or barcode that Guppy could have missed, as Guppy only removes adapters attached to the ends of the reads, not the ones in the middle of reads. 

`Script for Porechop`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  10
#SBATCH --job-name       Porechop
#SBATCH --mem            40G
#SBATCH --time           05:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load Porechop/0.2.4-gimkl-2020a-Python-3.8.2

porechop	-i /path_to_Nanolyse_output/Nanolyse/hectors_nanopore_merged_filtered.fastq.gz \ #Path to the output from Nanolyse
	-o /pat_for_Porechop_output/Porechop/hectors_nanopore_merged_filtered_porechop.fastq.gz \ #Path to the desired output for Porechop Results
	--threads ${SLURM_CPUS_PER_TASK} #Desired number of threads to use for the job, we will use the cpus-per-task assigned to the slurm script
```

## 5. Nanoqc
An extra quality assessment using Nanoqc to confirm the improvement in data quality after the previous procedures.
`Script for NanoQC`
```
#!/bin/bash -e

#SBATCH --nodes		1
#SBATCH --cpus-per-task	1
#SBATCH --ntasks		10
#SBATCH --job-name	NanoQC
#SBATCH --mem		10G
#SBATCH --time		05:00:00
#SBATCH --account	uoo02423
#SBATCH --output		%x.%j.out
#SBATCH --error		%x.%j.err
#SBATCH --hint		nomultithread

module purge
module load nanoQC/0.9.4-gimkl-2022a-Python-3.10.5

nanoQC	/path/to/Porecho/ouput/hectors_nanopore_merged_filtered_porechop.fastq.gz \
	-o /path/to/output/directory/NanoQCH2/HectorsNanoQC.html
```
## 6. Flye
We use Flye to assemble the genome with long reads. Flye assembly was performed with 3 iterations and with a 2.3GB genome size estimate as it is the estimate for a Dolphin species.

`Script for Flye`
```
#!/bin/bash -e

#SBATCH --nodes          1
#SBATCH --cpus-per-task  28
#SBATCH --job-name       Flye
#SBATCH --mem            900G
#SBATCH --time           5-00:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load Flye/2.9.1-gimkl-2022a-Python-3.10.5

flye --nano-hq /path_to_Porechop_output/hectors_nanopore_merged_filtered_porechop.fastq.gz \
     -o /path_for_Flye_assembly/FlyeHectors \
     -g 2.3G \ #Estimated genome size
     -t ${SLURM_CPUS_PER_TASK} \ #Number of threads as the assigned number os cpus-per-task for the slurm job
     -I 3 #Number of iterations
```
## 7. Quast
We use Quast to assess the genome assembly quality and compare the quality of the different genome assemblies (Hector's and Māui)

`Script for Quast`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  10
#SBATCH --job-name       Quast
#SBATCH --mem            20G
#SBATCH --time           05:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge 
module load QUAST/5.2.0-gimkl-2022a

quast.py /path_to_flye_assembly/flye_hectors_assembly.fasta \
	-t ${SLURM_CPUS_PER_TASK} \
	--eukaryote \ #Database against you would like to compare the BUSCO results of your assembly 
	--large \ #Size of your genome, large = typically > 100 Mbp)
	--conserved-genes-finding #BUSCO assessment of conserved genes in the lineage selected

-o /path_for_output/Quast
```
