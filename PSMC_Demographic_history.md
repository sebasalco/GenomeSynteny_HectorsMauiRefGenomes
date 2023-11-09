




#Files for PSMC

#PSMC bootstrap
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

psmc_plot.pl -R -u 2.3e-08 -g 12.5 -p hectors_bootstrap hectors_bootstrap.psmc
```


#Pseudo_PSMC
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
