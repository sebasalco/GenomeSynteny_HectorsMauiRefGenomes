# Genome Annotation
In the absence of RNAseq availability for the Hector’s and Māui dolphins, we annotated the genomes using GALBA and BRAKER3, with a reference-guided homology approach using the proteome and RNA-Seq of Tursiops truncatus. 
---
## 1. Repeatmasker
Prior to genome annotation, both genomes were screened for repetitive elements and softmasked.
`Script for Repeatmasker`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  16
#SBATCH --job-name       MaskHectors
#SBATCH --mem            75G
#SBATCH --time           24:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load RepeatMasker/4.1.0-gimkl-2020a

# Set the path to the Dfam HMM library
DFAM_LIB_PATH=/path/to/dfam/Dfam.hmm

RepeatMasker -species artiodactyl -xsmall -lib $DFAM_LIB_PATH hectors_genome.fasta

```
## 2. BRAKER3 annotation
Augustus config file
`Script for BRAKER3`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  24
#SBATCH --job-name       BrakerAnnotation
#SBATCH --mem            200G
#SBATCH --time           3-00:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load Singularity/3.11.3
module load AUGUSTUS/3.5.0-gimkl-2022a
module load BRAKER/3.0.3-gimkl-2022a-Perl-5.34.1

#export the path to augustus config copied above - prerequisites
export AUGUSTUS_CONFIG_PATH=/nesi/nobackup/uoo02423/Sebastian/AugConfigFile/config/

GENOME=/path/to/masked/genome/hectors_masked_genome.fasta
PROTEOME=/path/to/reference/proteome/Annotation/refdata/tursiops_proteins.fasta

srun braker.pl --threads=${SLURM_CPUS_PER_TASK} 
	--genome=$GENOME 
	--prot_seq=$PROTEOME


```
## 3. GALBA annotation.
`Script for GALBA`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  24
#SBATCH --job-name       GalbaAnnotation
#SBATCH --mem            100G
#SBATCH --time           7-00:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load Singularity/3.11.0
module load AUGUSTUS/3.5.0-gimkl-2022a

GENOME=/path/to/masked/genome/hectors_masked_genome.fasta
PROTEOME=/path/to/reference/proteome/Annotation/refdata/tursiops_proteins.fasta

export AUGUSTUS_CONFIG_PATH=/nesi/nobackup/uoo02423/Sebastian/AugConfigFile/config
export AUGUSTUS_BIN_PATH=/opt/nesi/CS400_centos7_bdw/AUGUSTUS/3.5.0-gimkl-2022a/bin

singularity run galba.sif galba.pl --version > galba.version
singularity exec --bind /path/to/working/directory/Annotation/GALBA/ galba.sif galba.pl \
	--species="Cephalorhynchus_hectori_maui" \
	--genome=$GENOME \
	--prot_seq=$PROTEOME \
	--AUGUSTUS_CONFIG_PATH=/nesi/nobackup/uoo02423/Sebastian/AugConfigFile/config/
	--threads ${SLURM_CPUS_PER_TASK} \
	--workingdir=/path/to/working/directory/Annotation/GALBA/ \
	--gff3 \
	--crf

```
## 4. TSEBRA merging annotations.
`Script for TSEBRA`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  24
#SBATCH --job-name       TsebraHec
#SBATCH --mem            100G
#SBATCH --time           4-00:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load TSEBRA/1.1.1-gimkl-2022a-Python-3.11.3

GALBA_GTF=/path/to/galba/output/gtf/Annotation/GALBA/GALBA/hectors_augustus.hints.gtf
GALBA_HINTS=/path/to/galba/output/hints_file/Annotation/GALBA/GALBA/hectors_galba_hintsfile.gff
BRAKER_GTF=/path/to/braker/output/gtf/Annotation/BRAKER/braker/hectors_braker.gtf
BRAKER_HINTS=/path/to/braker/output/hints_file/Annotation/Annotation/BRAKER/braker/hectors_hintsfile.gff

tsebra.py -g $GALBA_GTF,$BRAKER_GTF \
	-c default.cfg \
	-e $GALBA_HINTS,$BRAKER_HINTS \
	--filter_single_exon_genes \
	-o hec_combined_annot.gtf

```
## 5. Script to transform tsebra gtf file to proteins (doesn't need a slurm script as it is a low resources job)
```
module load AUGUSTUS/3.5.0-gimkl-2022a

gtf2aa.pl /path/to/masked/genome/hectors_masked_genome.fasta \
/path/to/tsebra/output/hec_combined_annot.gtf \
hectors_proteins_tsebra.fa

```
## 6. Gene Validator
`Script for GeneValidator`
```
#!/bin/bash -e


```
