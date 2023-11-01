# Merging nanopore and supernova assemblies and initial polishing steps

## 1. Merging the assemblies
We got two assemblies for each dolphin. In this step, we merged them to obtain an assembly for each dolphin and continue with the polishing steps separately for each dolphin.
---
## 1. Quickmerge
`Script for Quickmerge`
```
#!/bin/bash -e

#SBATCH --nodes		1
#SBATCH --cpus-per-task	1
#SBATCH --ntasks		10
#SBATCH --job-name	QuickmergeMF
#SBATCH --mem		100G
#SBATCH --time		48:00:00
#SBATCH --account	uoo02423
#SBATCH --output		%x_%j.out
#SBATCH --error		%x_%j.err
#SBATCH --hint		nomultithread

Install quick merge in conda environment

merge_wrapper.py  /nesi/nobackup/uoo02423/Sebastian/MergedAssemblies MauiFull\hectors_FlyeAssembly.fasta 
	         /nesi/nobackup/uoo02423/Sebastian/MergedAssemblies/MauiFull/hectors_supernova_assembly.fasta \
-hco 5.0 -c 1.5 -l 142457 -ml 5000 -p hectors_merged_assembly.fasta
```
## 2. Purgehaplotigs
We identified and resigned allelic contigs using purge_haplotigs to improve the genome assemblies.
`Script for purge_haplotigs`
```
#!/bin/bash -e

#SBATCH --cpus-per-task 	8
#SBATCH --job-name	purgehap
#SBATCH --partition	large,bigmem
#SBATCH --mem		50G
#SBATCH --time		10:00:00
#SBATCH --account	uoo02423
#SBATCH --output		%x_%j.out
#SBATCH --error		%x_%j.err
#SBATCH --hint		nomultithread

module purge
module load SAMtools/1.16.1-GCC-11.3.0
module load minimap2/2.24-GCC-11.3.0
module load BEDTools/2.30.0-GCC-11.3.0
module load purge_haplotigs/1.1.2-gimkl-2022a-Perl-5.34.1

# Commands below are run in sequences in the same working directory

# 1. Mapping longreads back to assembly
minimap2 -t 10 -ax map-ont \
	/path/to/merged/assembly/hectors_merged_assembly.fasta \
	/path/to/nanopore/long_reads/hectors_nanopore_filtered_reads.fastq.gz \
	--secondary=no | samtools sort -m 5G -o aligned.bam -T tmp.ali

# 2. Create histogram with purgehaplotigs
purge_haplotigs hist -t 10 -b aligned.bam \ # This bam file is produced in the step above
	-g /path/to/merged/assembly/hectors_merged_assembly.fasta

# 3. Coverage step with custom parameters
purge_haplotigs cov -i aligned.bam.gencov -l 1 -m 11 -h 190 -o coverage_stats.csv

# 3. Purging step
purge_haplotigs purge -g /path/to/merged/assembly/hectors_merged_assembly.fasta \
	-c coverage_stats.csv
	-b aligned.bam -d
```
## 3. Rails & Cobbler.
Rails and cobbler was used as the first gap-filling step as well as an initial scaffolding step
`Script for mkoutput`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  12
#SBATCH --job-name       RailsCobbler
#SBATCH --mem            120G
#SBATCH --time           2-00:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load Perl/5.34.1-GCC-11.3.0
module load minimap2/2.24-GCC-11.3.0
module load SAMtools/1.13-GCC-9.2.0

sh runRAILSminimapSTREAM.sh hectors_curated.fasta \ #Output file from purge_haplotigs
	/path/to/nanopore/long_reads/hectors_nanopore_filtered_reads.fastq.gz 250 0.85 500 2 ont
	/path/to/samtools/opt/nesi/CS400_centos7_bdw/SAMtools/1.13-GCC-9.2.0/bin/samtools
```
## 4. LR_GapCloser
We finished the gap closing using LR_Gapcloser and the nanopore long reads, performing a total of 10 gap-closing iterations
`Script for LR_Gapcloser`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  12
#SBATCH --job-name       lr_gapc
#SBATCH --mem            100G
#SBATCH --time           24:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load BWA/0.7.17-GCC-11.3.0
export PATH=/nesi/nobackup/uoo02423/Sebastian/MergedAssemblies/Hectors42x/LR_Gapcloser/src/:$PATH

sh LR_Gapcloser.sh -i /path/to/output/from/rails&cobbler/hectors_RandC_output.fasta \ #Output file from purge_haplotigs
	-l /path/to/nanopore/long_reads/hectors_nanopore_filtered_reads.fastq.gz \
	-s n 
	-t ${SLURM_CPUS_PER_TASK}
	-r 10
```
