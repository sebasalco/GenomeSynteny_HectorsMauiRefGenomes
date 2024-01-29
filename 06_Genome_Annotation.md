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
Gene validator is a program with a high-resource usage. The program is installed easily cloning the repository from github (https://github.com/wurmlab/genevalidator). For the gene validation we constructed a diamond database '4cetcombined.dmnd' of the annotated genes for 4 cetaceans (Orca, bottlenose, vaquita and beluga). Then, the gene validation was performed against that constructed database. 
`Script for GeneValidator`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  64
#SBATCH --job-name       gv_hectors
#SBATCH --mem            550G
#SBATCH --time           24:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load DIAMOND/2.1.6-GCC-11.3.0

export PATH=/nesi/nobackup/uoo02423/Sebastian/MergedAssemblies/Hectors42x/Annotation/genevalidator/:$PATH

diamond makedb --db /nesi/nobackup/uoo02423/Sebastian/MergedAssemblies/Hectors42x/Annotation/4cetcombined.dmnd \
    --in /nesi/nobackup/uoo02423/Sebastian/MergedAssemblies/Hectors42x/Annotation/refdata/4cet_combined_proteins.fasta

diamond blastp --db 4cetcombined.dmnd -q /nesi/nobackup/uoo02423/Sebastian/MergedAssemblies/Hectors42x/Annotation/2hectors_proteins_tsebra.fa \
    --more-sensitive --outfmt 5 --threads ${SLURM_CPUS_PER_TASK} > 2hectors_proteins.diamond.xml

genevalidator -o 2hectors_gv -x 2hectors_proteins.diamond.xml \
    --raw_sequences /nesi/nobackup/uoo02423/Sebastian/MergedAssemblies/Hectors42x/Annotation/refdata/4cet_combined_proteins.fasta -n ${SLURM_CPUS_PER_TASK} \
    -m ${SLURM_CPUS_PER_TASK} /nesi/nobackup/uoo02423/Sebastian/MergedAssemblies/Hectors42x/Annotation/hectors_proteins_tsebra.fa
```
## 7. eggNOG mapper 
Previous to this step all the protein and gtf files had a generic ID for each entry (e.g. anno2.g29559.t1). We used eggnog-mapper to asing the gene IDS to each of the transcripts and proteins. Eggnog output is a tab delimited table with all of our generic IDs for protein and their new gen ID. Eggnog has its own database that needs to be downloaded previous to the analysis.
`Script for eggnog-mapper`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  24
#SBATCH --job-name       eggnog
#SBATCH --mem            50G
#SBATCH --time           24:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load eggnog-mapper/2.0.1b-gimkl-2020a-Python-2.7.18

DATA_PATH=/opt/nesi/db/eggnog_db/data/

echo $DATA_PATH

emapper.py -m diamond \
     -i /path/to/transformed/tsebra/output/hectors_proteins_tsebra.fa \
     -o hectors_annotated --data_dir $DATA_PATH
```
To assign the new gene IDs from the tab delimited output to the protein file
```
#!/bin/bash -e

#SBATCH --cpus-per-task  4
#SBATCH --job-name       Id_reassign
#SBATCH --mem            10G
#SBATCH --time           2:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

# Create an associative array to store gene information
declare -A gene_info

# Read gene information from the text file
while read -r gene_id gene_name; do
    gene_info["$gene_id"]=$gene_name
done < hectors_final_annotIDfile.txt

# Process the protein FASTA file
while IFS= read -r line; do
    if [[ $line == '>'* ]]; then
        protein_id=${line#>}

        # Look up gene name from the associative array
        gene_name=${gene_info["$protein_id"]}

        # If gene name exists, update the header; otherwise, use the original header
        if [ -n "$gene_name" ]; then
            echo ">$protein_id $gene_name"
        else
            echo "$line"
        fi
    else
        echo "$line"
    fi
done < hectors_proteins_tsebra.fa > final_hectors_proteins_tsebra.fa
```
