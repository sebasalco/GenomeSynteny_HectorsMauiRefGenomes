# Genome Wide Heterozygosity Analysis using angsd (low coverage) and VCFtools
---
## 1. Calculating GWH using angsd. 
First we mapped the WGS short reads of each dolphin to their own reference genome, obtaining the bam files. Same procedure as the one performed in step 8 for the PSMC.

`Script for angsd GWH`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  12
#SBATCH --job-name       TestAngsd
#SBATCH --mem            150G
#SBATCH --time           5:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load angsd/0.935-GCC-9.2.0

GENOME=/path/to/genome/hectors_aut_genome.fasta
BAM=/path/to/bam/file/hectors.aligned.sorted.bam

angsd -i $BAM -anc $GENOME -dosaf 1 -GL 1

realSFS angsdput.saf.idx > est.ml

```

## 1. Calculating GWH in sliding windows. 
We extracted the mean GWH per chromosome and we calculated GWH in 10kb sliding windows across the whole genome filtering sites with low quality, low depth, and removing indels.

`Script for 10kb sliding window GWH`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  12
#SBATCH --job-name       winAngsd
#SBATCH --mem            200G
#SBATCH --time           24:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load angsd/0.935-GCC-9.2.0

GENOME1='/nesi/nobackup/uoo02423/Sebastian/HectorsMauiGenomeFiles/hectors_aut_genome.fasta'
GENOME2='/nesi/nobackup/uoo02423/Sebastian/HectorsMauiGenomeFiles/maui_aut_genome.fasta'
BAM1='hectors.bam'
BAM2='maui.bam'

angsd -i hectors.aligned.sorted.bam -anc $GENOME1 -ref $GENOME1 -dosaf 1 \
-gl 2 -out hectors -nthreads ${SLURM_CPUS_PER_TASK} -minQ 20

angsd -i maui.aligned.sorted.bam -anc $GENOME2 -ref $GENOME2 -dosaf 1 \
-gl 2 -out maui -nthreads ${SLURM_CPUS_PER_TASK} -minQ 20

realSFS hectors.saf.idx

realSFS maui.saf.idx

thetaStat do_stat hectors.saf.idx -win 50000 -step 5000 -outnames hec_theta.thetasWindow.gz

thetaStat do_stat maui.saf.gz -win 50000 -step 5000 -outnames maui_theta.thetasWindow.gz
```




