# Historical demography of the Hector’s and Māui dolphin constructed using the Pairwise Sequentially Markovian Coalescent (PSMC) model
---
## 1. Files for PSMC. 
First we mapped the WGS short reads of each dolphin to their own reference genome.

`Script for minimap2`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  8
#SBATCH --job-name       Map_2Genome
#SBATCH --mem            25G
#SBATCH --time           16:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module load SAMtools/1.16.1-GCC-11.3.0
module load minimap2/2.24-GCC-11.3.0

GENOME=/nesi/nobackup/uoo02423/Sebastian/MergedAssemblies/Hectors42x/CHROMOSOME_GENOME/hectors_aut_genome.fasta

minimap2 -ax sr -t ${SLURM_CPUS_PER_TASK} \
$GENOME \
/path/to/WGSfiles/WGS-Che11CB067_combined_trimmed_R1.fastq.gz \
/path/to/WGSfiles/WGS-Che11CB067_combined_trimmed_R2.fastq.gz \
| samtools sort -@10 -O BAM -o wgs_hectors_Assembly.bam$
```

We indexed the genome, and used the output .bam file to call the variants with bcftools mpileup

`Script for bcftools mpileup`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  12
#SBATCH --job-name       FilexPSMCHec
#SBATCH --mem            15G
#SBATCH --time           24:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load psmc/0.6.5-gimkl-2018b
module load SAMtools/1.16.1-GCC-11.3.0
module load BCFtools/1.16-GCC-11.3.0

GENOME=/nesi/nobackup/uoo02423/Sebastian/MergedAssemblies/Hectors42x/CHROMOSOME_GENOME/hectors_aut_genome.fasta

samtools faidx $GENOME

samtools index wgs_hectors_Assembly.bam

bcftools mpileup -Q 30 -q 30 -O v \
-f $GENOME wgs_hectors_Assembly.bam |  bcftools call -c | vcfutils.pl vcf2fq -d 5 -D 34 -Q 30 > hectors_forpsmc.fq
```

## 2. PSMC bootstrap
Here we transforme de .fq file to psmc and performed a bootstrap of 100.

`Script for PSMC bootstrap`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  6
#SBATCH --job-name       Bt_PSMCHec
#SBATCH --mem            10G
#SBATCH --time           01:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load psmc/0.6.5-gimkl-2018b

fq2psmcfa hectors_forpsmc.fq > hectors.psmcfa

splitfa hectors.psmcfa > split.psmcfa

psmc -N25 -t15 -r5 -d -p "8+23*2+9+1" -o hectors.psmc hectors.psmcfa

seq 100 | xargs -i echo psmc -N25 -t15 -r5 -b -p "8+23*2+9+1" -o round-{}.psmc split2m.psmcfa | sh

cat hectors2m.psmc round-*.psmc > hectors_bootstrap.psmc

psmc_plot.pl -R -u 1.3e-08 -g 12.5 -p hectors_bootstrap hectors_bootstrap.psmc
```


## 3. Pseudo_PSMC
To analyze when Hector’s and Māui dolphins could have been reproductively isolated (cessation of gene flow) from each other, we performed a pseudodiploid/hybrid PSMC using seqtk.

`Script for pseudo PSMC`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  6
#SBATCH --job-name       Pseudo_PSMC
#SBATCH --mem            15G
#SBATCH --time           12:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load psmc/0.6.5-gimkl-2018b
module load seqtk/1.4-GCC-11.3.0

#Generate a merged psmcfa file
 
seqtk mergefa -rq20 hectors_forpsmc.fq maui_forpsmc.fq | fq2psmcfa -q30 - > merged.psmcfa

#Generate psmc file
 
psmc -N25 -t15 -r5 -p "8+23*2+9+1" -o merged.psmc merged.psmcfa
 
#Generate plot
 
smc_plot.pl -R -u 1.5e-08 -g 12.5 -p 1.5merged_plot merged.psmc
```
