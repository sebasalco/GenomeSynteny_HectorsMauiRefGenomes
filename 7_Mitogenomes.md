# Extracting mitogenomes from short-read WGS data using MitoFinder.
---
## 1. MitoFinder
We use MitoFinder with the megahit assembler option and using the Commerson's dolphin mitogenome as reference. 

`Script for MitoFinder`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  12
#SBATCH --job-name       Hec_mito
#SBATCH --mem            50G
#SBATCH --time           48:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module load BLAST/2.13.0-GCC-11.3.0
module load MEGAHIT/1.2.9-gimkl-2022a-Python-3.10.5
module load Singularity/3.11.3

mitofinder_v1.4.1.sif --megahit -j Hectors_mitogenome \
-1 WGS-Che11CB067_combined_trimmed_R1.fastq.gz \
-2 WGS-Che11CB067_combined_trimmed_R2.fastq.gz \
-r Comm_mito.gb \
-o 2 \
```

## 2. Assesing coverage 
First we mapped the reads to the obtained mitogenomes using Minimap2

`Script for Minimap`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  8
#SBATCH --job-name       C_MitoCovMaui
#SBATCH --mem            25G
#SBATCH --time           16:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module load SAMtools/1.16.1-GCC-11.3.0
module load minimap2/2.24-GCC-11.3.0

minimap2 -ax sr -t ${SLURM_CPUS_PER_TASK} \
/nesi/nobackup/uoo02423/Sebastian/HectorsMauiGenomeFiles/Mitogenomes/Hectors_mitogenome.fasta \
/nesi/nobackup/uoo02423/Sebastian/HectorsMauiGenomeFiles/Mitogenomes/WGS-Che11CB067_combined_trimmed_R1.fastq.gz \
/nesi/nobackup/uoo02423/Sebastian/HectorsMauiGenomeFiles/Mitogenomes/WGS-Che11CB067_combined_trimmed_R2.fastq.gz \
| samtools sort -@10 -O BAM -o mito_hectors.bam$
```
We indexed the bam files and used bedtools to obtain a bedgraph with the coverage

`Script for bedtools`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  8
#SBATCH --job-name       bedCoverage
#SBATCH --mem            25G
#SBATCH --time           12:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module load BEDTools/2.30.0-GCC-11.3.0
module load SAMtools/1.16.1-GCC-11.3.0

samtools index mito_hectors.bam

bedtools genomecov -ibam mito_hectors.bam -bg -g Hectors_mitogenome.bed > hectors_coverage.bedgraph
```
