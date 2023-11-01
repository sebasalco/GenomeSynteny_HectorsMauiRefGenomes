# Second polishing stage
A second polishing stage was performed to the chromosome-level genomes obtained from the previous scaffold stage
---
## 1. Purgehaplotigs
`Script for purgehaplotigs`
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

# 1. Mapping longreads back to the chromosome-level assemblies
minimap2 -t 10 -ax map-ont \
	/path/ragtag/output/hectors_genome.fasta \
	/path/to/nanopore/long_reads/hectors_nanopore_filtered_reads.fastq.gz \
	--secondary=no | samtools sort -m 5G -o aligned.bam -T tmp.ali

# 2. Create histogram with purgehaplotigs
purge_haplotigs hist -t 10 -b aligned.bam \ # This bam file is produced in the step above
	-g /path/ragtag/output/hectors_genome.fasta

# 3. Coverage step with custom parameters
purge_haplotigs cov -i aligned.bam.gencov -l 1 -m 11 -h 190 -o coverage_stats.csv

# 3. Purging step
purge_haplotigs purge -g /path/ragtag/output/hectors_genome.fasta \
	-c coverage_stats.csv
	-b aligned.bam -d
```
## 2. BLAST for BlobTools
We blasted the .fasta file output from the previous purge_haplotigs stage agains the BLAST database.
`Script for BLAST`
```
#!/bin/bash -e

#SBATCH --partition      milan
#SBATCH --cpus-per-task  64
#SBATCH --job-name       BlastHectors
#SBATCH --mem            800G
#SBATCH --time           48:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load BLAST/2.13.0-GCC-11.3.0
module load BLASTDB/2023-01

# This script takes one argument, the FASTA file of query sequences.
QUERIES=/path/to/output/from/purging/hectors_genome_curated.fasta
FORMAT="6 qseqid staxids bitscore std evalue"
BLASTOPTS="-task blastn -evalue 1e-25 -max_target_seqs 10 -max_hsps 1"
BLASTAPP=blastn
DB=nt

# Keep the database in RAM
cp $BLASTDB/{$DB,taxdb}* $TMPDIR/
export BLASTDB=$TMPDIR

$BLASTAPP $BLASTOPTS -db $DB -query $QUERIES -outfmt "$FORMAT" \
-out Hectors_blastn.out -num_threads $SLURM_CPUS_PER_TASK
```
## 3. Diamond for BlobTools.
`Script for Diamond`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  4
#SBATCH --job-name       Diamond
#SBATCH --mem            50G
#SBATCH --time           48:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load DIAMOND/2.0.6-GCC-9.2.0

diamond blastx \
              --query /path/to/output/from/purging/hectors_genome_curated.fasta \
              --db path/to/the/downloaded/and/formatted/uniprot/database/reference_proteomes.dmnd \
              --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
              --sensitive \
              --max-target-seqs 1 \
              --evalue 1e-25 \
              --threads $SLURM_CPUS_PER_TASK \
              --out Hectors_diamond_blastx.out
```
## 4. BlobTools.
`Script for BlobTools`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  12
#SBATCH --job-name       BlobHectors
#SBATCH --mem            40G
#SBATCH --time           07:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

blobtools create \
	--fasta /path/to/output/from/purging/hectors_genome_curated.fasta \
	--meta Hectors_database.yaml --taxid 1183062 \
	--taxdump /path/to/database/Database/taxdump \
	${PWD}/DATASETS/Hectors_assembly

blobtools add \
	--hits path/to/blastn/hits/Hectors_blastn.out \
	--hits path/to/diamond.blastx/hits/Hectors_diamond_blastx.out \
	--busco /path/to/cetartiodactyla/database/run_cetartiodactyla_odb10/full_table.tsv \
	--cov /path/to/coverage/Hectors_Assembly.bam \
	${PWD}/DATASETS/Hectors_assembly


blobtools add \
	--hits path/to/blastn/hits/EW_Assembly.ncbi.blastn.out \
	--hits path/to/diamond.blastx/hits/EW_Assembly.diamond.blastx.out \
	--cov path/to//coverage/bam/file/EW_Assembly.bam \
	--taxrule bestsumorder \
	--taxdump /path/to/database/Database/taxdump \
	${PWD}/DATASETS/Hectors_assembly

blobtools filter \
        --query-string "length--Min=1000&Hectors_Assembly.bam--Min=5.00" \
        --fasta /path/to/output/from/purging/hectors_genome_curated.fasta \
        --output path/to/output/folder/filter \
        ${PWD}/DATASETS/Hectors_assembly
```
## 5. Pilon.
Error correction was performed using Pilon with short-read sequences.
`Script for Pilon`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  16
#SBATCH --job-name       Pilon
#SBATCH --mem            50G
#SBATCH --time           24:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load Pilon/1.24-Java-15.0.2

java -Xmx16G -jar $EBROOTPILON/pilon.jar \
    --genome /path/to/output/from/purging/hectors_genome_curated.fasta \
    --frags illumina_trimmed_mapped_hectors_genome.sorted.bam \
    --changes --vcf --diploid --threads ${SLURM_CPUS_PER_TASK} \
    --output hectors_pilon_polish1 \
    --minmq 30 \
    --minqual 30

```
## 6. Final Scaffolding
The final orientating of chromosomes was performed against the bottlenose dolphin as it is the closest species chromosome-level genome. This final scaffolding was performed using NUCMER.
`Script for Ragtag`
```
#!/bin/bash -e

#SBATCH --cpus-per-task 	8
#SBATCH --job-name	RagtagScff
#SBATCH --mem		50G
#SBATCH --time		10:00:00
#SBATCH --account	uoo02423
#SBATCH --output		%x_%j.out
#SBATCH --error		%x_%j.err
#SBATCH --hint		nomultithread

module purge
module load minimap2/2.24-GCC-11.3.0

export PATH=/path/to/ragtag/installation/RagTag/:$PATH

ragtag.py scaffold /path/to/reference/RefCetaceanGenomes/Tursiops/T_truncatus.fasta \
/path/output/from/pilon/hectors_pilon_polish.fasta \
-aligner 'nucmer'
```
