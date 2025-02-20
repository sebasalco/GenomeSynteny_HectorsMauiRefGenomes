# We use FCSGX screening to match with NCBI standars on contamination and adaptor removal. For this, we downloaded the FCSGX database and the singularity image of the program. We used a temporary directory in RAM to store the dataset for efficiency.

## 1. FCSGX screen and clean
`Script for FCSGX`
```
#!/bin/bash -e

#SBATCH --cpus-per-task 12
#SBATCH --job-name	scrfgcgx
#SBATCH --mem		500GB
#SBATCH --time		12:00:00
#SBATCH --account	uoo02423
#SBATCH --output	%x_%j.out
#SBATCH --error		%x_%j.err
#SBATCH --hint		nomultithread

module purge
module load Singularity/3.11.3

export FCS_DEFAULT_IMAGE=fcs-gx.sif

# Set the database
DB_ORIG=/nesi/nobackup/uoo02423/Sebastian/SyntenyRev/FCSGX/databaseFCSGX
DB_RAM=$TMPDIR/databaseFCSGX

# Check available TMPDIR space before copying
DB_SIZE=$(du -sh $DB_ORIG | cut -f1)
TMP_SIZE=$(df -h $TMPDIR | tail -1 | awk '{print $4}')
echo "Database size: $DB_SIZE | TMPDIR available: $TMP_SIZE"

cp -r $DB_ORIG $DB_RAM  # Copy database to RAM
export GXDB_LOC=$DB_RAM  # Set database location for FGS-GX

# Set query genome file & output
QUERIES=/nesi/nobackup/uoo02423/Sebastian/SyntenyRev/MerQV/Hectors_MergedAssembly.fasta.gz  # Ensure it's gzipped
OUTDIR=/nesi/nobackup/uoo02423/Sebastian/SyntenyRev/FCSGX/output
mkdir -p $OUTDIR

# Run Genome Screening inside Singularity
python3 ./fcs.py screen genome --fasta $QUERIES --out-dir $OUTDIR --gx-db "$GXDB_LOC/gxdb" --tax-id 37035

python3 ./fcs.py clean genome --action-report hectors_contamination_report.txt --output hectors_clean.fasta --contam-fasta-out contam.fasta
```
