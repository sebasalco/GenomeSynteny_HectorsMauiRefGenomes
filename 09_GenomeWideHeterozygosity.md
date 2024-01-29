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

#SBATCH --cpus-per-task  8
#SBATCH --job-name       HecHeter
#SBATCH --mem            20G
#SBATCH --time           5:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load SAMtools/1.16.1-GCC-11.3.0
#module load BCFtools/1.16-GCC-11.3.0
module load VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1

bcftools query -f '%CHROM\t[%GT]\n' hectors_final_variants.vcf | awk 'BEGIN{OFS="\t"} {chrom = substr($1, 10); if ($2 ~ /[0-9]\/[0-9]/) het[chrom]++} END{for (chrom in het) print chrom, het[chrom]/NR}' > heterozygosity_per_chromosome.txt

vcftools --vcf hectors_final_variants.vcf --window-pi 10000 --out filter_hectors_10kb \
  --minQ 30 \
  --minDP 3 \
  --maxDP 1000 \
  --max-missing 0.2 \
  --maf 0.05 \
  --remove-indels
```




