# We evaluated reads mapping quality before and after the reference scaffolding to assess if the scaffolding approach has any effect in breaking long-read connections. We aligned raw ONT reads to the prescaffolded genomes and to the final genome assemblies.

## 1. Minimap2 indexing 
`Script for minimap2`
```
#!/bin/bash -e

#SBATCH --cpus-per-task 12
#SBATCH --job-name      mmp2ind
#SBATCH --mem           50G
#SBATCH --time          06:00:00
#SBATCH --account       uoo02423
#SBATCH --output        %x_%j.out
#SBATCH --error         %x_%j.err
#SBATCH --hint          nomultithread

module purge
module load minimap2/2.28-GCC-12.3.0

minimap2 -d Hectors_merged_Assembly.mmi Hectors_merged_Assembly.fasta #Pre scaffolded genome assembly

minimap2 -d merged_MauiFull_Assembly.mmi merged_MauiFull_Assembly.fasta #Pre scaffolded genome assembly

minimap2 -d hectors_ordered.mmi hectors_genome.fasta #Reference scaffolded genome assembly

minimap2 -d maui_ordered.mmi maui_genome.fasta #Reference scaffolded genome assembly
```
## 2. Minimap2 alignment 
`Script for minimap2`
```
#!/bin/bash -e

#SBATCH --cpus-per-task 12
#SBATCH --job-name      alnONT
#SBATCH --mem           50G
#SBATCH --time          06:00:00
#SBATCH --account       uoo02423
#SBATCH --output        %x_%j.out
#SBATCH --error         %x_%j.err
#SBATCH --hint          nomultithread

module purge
module load minimap2/2.28-GCC-12.3.0
module load SAMtools/1.19-GCC-12.3.0

minimap2 -ax map-ont merged_MauiFull_Assembly.mmi Maui6runsFiltered.fastq.gz | samtools sort \ #Map ONT reads to pre-scaffolded assembly 
-o ont_maui_non_scaffolded.bam
samtools index ont_maui_non_scaffolded.bam

minimap2 -ax map-ont maui_genome.mmi Maui6runsFiltered.fastq.gz | samtools sort \ #Map ONT reads to referece scaffolded final assembly 
-o ont_maui_scaffolded.bam
samtools index ont_maui_scaffolded.bam

minimap2 -ax map-ont Hectors_merged_Assembly.mmi Hector2runsFiltered.fastq.gz.fastq.gz | samtools sort \ #Map ONT reads to pre-scaffolded assembly 
-o ont_hec_non_scaffolded.bam
samtools index ont_hec_non_scaffolded.bam

minimap2 -ax map-ont hectors_genome.mmi Hector2runsFiltered.fastq.gz.fastq.gz | samtools sort \ #Map ONT reads to referece scaffolded final assembly 
-o ont_hec_scaffolded.bam
samtools index ont_hec_scaffolded.bam
```
## 3. Samtools flagstat alignment evaluation 
`Script for samtools`
```
#!/bin/bash -e

#SBATCH --cpus-per-task 12
#SBATCH --job-name      stats
#SBATCH --mem           10G
#SBATCH --time          06:00:00
#SBATCH --account       uoo02423
#SBATCH --output        %x_%j.out
#SBATCH --error         %x_%j.err
#SBATCH --hint          nomultithread

module purge
module load SAMtools/1.19-GCC-12.3.0

samtools flagstat ont_hec_scaffolded.bam
samtools flagstat ont_hec_non_scaffolded.bam
samtools flagstat ont_maui_scaffolded.bam
samtools flagstat ont_maui_non_scaffolded.bam

```
## 4. We evaluated the completeness and quality of the final assemblies using the kmer based statistics from Merqury QV 
`Script for meryl, create database from WGS`
```
#!/bin/bash -e

#SBATCH --cpus-per-task 12
#SBATCH --job-name      merylkd
#SBATCH --mem           50G
#SBATCH --time          06:00:00
#SBATCH --account       uoo02423
#SBATCH --output        %x_%j.out
#SBATCH --error         %x_%j.err
#SBATCH --hint          nomultithread

module purge
module load Merqury/1.3-Miniconda3

meryl k=21 count WGS-Che11CB067_combined_trimmed_R1.fastq.gz output kmerdb_hec_wgsr1.meryl

meryl k=21 count WGS-Che11CB067_combined_trimmed_R2.fastq.gz output kmerdb_hec_wgsr2.meryl

meryl union-sum output hec.meryl kmerdb_hec_wgsr1.meryl kmerdb_hec_wgsr2.meryl
```
`Script for merqury`
```
#!/bin/bash -e

#SBATCH --cpus-per-task 12
#SBATCH --job-name      qvmerq
#SBATCH --mem           50G
#SBATCH --time          06:00:00
#SBATCH --account       uoo02423
#SBATCH --output        %x_%j.out
#SBATCH --error         %x_%j.err
#SBATCH --hint          nomultithread

# Load required modules
module purge
module load Merqury/1.3-Miniconda3
module load BEDTools/2.31.1-GCC-12.3.0
module load SAMtools/1.19-GCC-12.3.0

# Ensure MERQURY environment variable is set
export MERQURY=/opt/nesi/CS400_centos7_bdw/Merqury/1.3-Miniconda3/merqury

# Run Merqury
$MERQURY/merqury.sh hec.meryl \
/path/to/genome/hectors_genome.fasta qv_hec_fin
```
