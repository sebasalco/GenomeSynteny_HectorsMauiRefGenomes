# 10x Chromium linked-reads assembly. 

## 1. Extracting the data
We got `*.bcl` data back from the sequencer, so first needed to generate the fastq files based off these. We did this separately for each lane after first unzipping (2hrs and 3GB of RAM was enough):
```
tar -zvxf 200702_A00488_0073_BHTL7YDMXX-lane1.tar.gz
tar -zvxf 200702_A00488_0073_BHTL7YDMXX-lane2.tar.gz
```

## 2. Running the supernova assembly
We ran our assemblies a few different ways, an initial 'maui' and 'hectors' run to gauge the depth effective run depth etc. Based on these, additional runs were done at 56x raw depth, and 42x effective depth. Following the initial supernova run, we extracted the two pseudohaploid genomes for each run using mkoutput files.
---
## 1. Supernova makefastq
`Script for makefastq`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  20
#SBATCH --job-name       makefastq
#SBATCH --mem            105G
#SBATCH --time           10:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load bcl2fastq2/2.20.0-gimkl-2018b
module load Supernova/2.1.1

supernova mkfastq --run=200702_A00488_0073_BHTL7YDMXX-lane1 --id=lane1_fastq --samplesheet=10x_samplesheet.csv --qc --jobmode=local --localcores=36 --localmem=105 --ignore-dual-index 

supernova mkfastq --run=200702_A00488_0073_BHTL7YDMXX-lane2 --id=lane2_fastq --samplesheet=10x_samplesheet_L2.csv --qc --jobmode=local --localcores=36 --localmem=105 --ignore-dual-index 
```

## 2. Supernova makeoutput
`Script for mkoutput`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  12
#SBATCH --job-name       PycoQC
#SBATCH --mem            400G
#SBATCH --time           7-00:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load pycoQC/2.5.2-gimkl-2020a-Python-3.8.2

module load Supernova/2.1.1

supernova run --id=H_10xSN \
              --fastqs=/path/to/linked-read/fastq/files

supernova mkoutput \
        --style=pseudohap \
        --asmdir=path/to/working-directory-of-supernova/H_10xSN/outs/assembly \ #Path to the output directory 'assembly' created by supernova
        --outprefix=./hectors_supernova_assembly.fasta #Prefix for output
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

quast.py /path_to_supernova_assembly/hectors_supernova_assembly.fasta \
	-t ${SLURM_CPUS_PER_TASK} \
	--eukaryote \ #Database against you would like to compare the BUSCO results of your assembly 
	--large \ #Size of your genome, large = typically > 100 Mbp)
	--conserved-genes-finding #BUSCO assessment of conserved genes in the lineage selected

-o /path_for_output/Quast
```
