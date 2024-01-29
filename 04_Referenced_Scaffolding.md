# Reference scaffolding
Reference scaffolding against the reference genomes published for Commerson's dolphin, Vaquita, Bottlenose dolphin, Orca and Blue Whale. Each step was performed for Hector's and Māui dolphin assemblies separately 

## 1. Ragtag Scaffold
The output from LR_Gapcloser was used as query genome for this scaffolding step, a scaffolding step was performed using that same input with the reference genomes separately obtaining the gtf output files for the final scaffolding with the info of all these 5 separate scaffoldings.
---
## 1. Ragtag Scaffold
`Script for Ragtag Scaffold`
```
#!/bin/bash -e

#SBATCH --cpus-per-task 8
#SBATCH --job-name	RagtagScff
#SBATCH --mem		50G
#SBATCH --time		10:00:00
#SBATCH --account	uoo02423
#SBATCH --output	%x_%j.out
#SBATCH --error		%x_%j.err
#SBATCH --hint		nomultithread

module purge
module load minimap2/2.24-GCC-11.3.0

export PATH=/path/to/ragtag/installation/RagTag/:$PATH

ragtag.py scaffold /path/to/reference/RefCetaceanGenomes/Commersonii/Commersonii.fasta \
/path/output/from/LR_Gapcloser/hectors_gapclosed_scf.fasta

ragtag.py scaffold /path/to/reference/RefCetaceanGenomes/Tursiops/T_truncatus.fasta \
/path/output/from/LR_Gapcloser/hectors_gapclosed_scf.fasta

ragtag.py scaffold /path/to/reference/RefCetaceanGenomes/Vaquita/P_sinus.fasta \
/path/to/assembly/output/from/LR_Gapcloser/hectors_gapclosed_scf.fasta

ragtag.py scaffold /path/to/reference/RefCetaceanGenomes/Orca/O_orca.fasta \
path/to/assembly/output/from/LR_Gapcloser/hectors_gapclosed_scf.fasta

ragtag.py scaffold /path/to/reference/RefCetaceanGenomes/BlueWhale/B_musculus.fasta \
path/to/assembly/output/from/LR_Gapcloser/hectors_gapclosed_scf.fasta
```
## 2. Merging Scaffolds
In this step, we merged the results from the 5 different scaffolding performed in the previous step. Agp files encode adjacency and gap information and can be hierarchized by importance (weight). In this case, we did not assign with to the agp files.
`Script for Ragtag merge`
```
#!/bin/bash -e

#SBATCH --cpus-per-task 8
#SBATCH --job-name	RagtagMerge
#SBATCH --mem		50G
#SBATCH --time		10:00:00
#SBATCH --account	uoo02423
#SBATCH --output	%x_%j.out
#SBATCH --error		%x_%j.err
#SBATCH --hint		nomultithread

module purge
module load minimap2/2.24-GCC-11.3.0

export PATH=/path/to/ragtag/installation/RagTag/:$PATH 

ragtag.py merge path/to/assembly/output/from/LR_Gapcloser/hectors_gapclosed_scf.fasta \ #Query genome, output from LR_gapcloser
	hectors_scaffold_ref_Oorca.agp \ #Agp file resulted from previous scaffolding against the Orca
	hectors_scaffold_ref_Psinus.agp \ #Agp file resulted from previous scaffolding against the Vaquita
	hectors_scaffold_ref_Ttruncatus.agp \ #Agp file resulted from previous scaffolding against the Bottlenose
	hectors_scaffold_ref_Ccommersonni.agp \ #Agp file resulted from previous scaffolding against the Commerson's
	hectors_scaffold_ref_Bmusculus.agp #Agp file resulted from previous scaffolding against the Blue whale
```
## 3. SANS phylogenomic tree.
We constructed a phylogenetic tree to confirm the phylogenetic relationships of the genomes obtained. This tree was constructed using the output chromosome-level assemblies obtained using Ragtag and the reference genomes used for that scaffolding.
`Script for SANS`
```
#!/bin/bash -e

#SBATCH --cpus-per-task	16
#SBATCH --job-name	SANSCet
#SBATCH --mem		300G
#SBATCH --time		72:00:00
#SBATCH --account	uoo02423
#SBATCH --output	%x_%j.out
#SBATCH --error		%x_%j.err
#SBATCH --hint		nomultithread

export PATH=/path/to/sans/installation/directory/sans/:$PATH

SANS -i CetaceanList.txt \ #List of the cetacean reference genomes paths
     -N CetaceanTree.newick 
     -T ${SLURM_CPUS_PER_TASK} 
     -f tree #Output type
     -b 1000 #Bootstrap
     -C # Consensus tree
```
`CetaceanList.txt`
```
/path/to/reference/RefCetaceanGenomes/Commersonii/Cephalorhynchus_commersonii.fasta
/path/to/reference/RefCetaceanGenomes/Orca/Orcinus_orca.fasta
/path/to/reference/RefCetaceanGenomes/Vaquita/Phocoena_sinus.fasta
/path/to/reference/RefCetaceanGenomes/Tursiops/Tursiops_truncatus.fasta
/path/to/reference/RefCetaceanGenomes/BlueWhale/Balaenoptera_musculus.fasta
/path/ragtag/output/hectors_genome.fasta
/path/ragtag/output/maui_genome.fasta
```
## 4. miniBUSCO assessment.
We performed a miniBUSCO analysis of the chromosome-level genomes obtained at the end of this stage.
`Script for miniBUSCO`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  24
#SBATCH --job-name       minibuscoHec
#SBATCH --mem            50G
#SBATCH --time           12:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load miniBUSCO/0.2.1-gimkl-2022a
module load Python/3.11.3-gimkl-2022a

minibusco.py run -a /path/ragtag/output/hectors_genome.fasta \
-t ${SLURM_CPUS_PER_TASK} \
-o H_minibusco \
-l cetartiodactyla #Lineage \
-L /path/to/busco/lineage/directory/Busco/Lineage
```
