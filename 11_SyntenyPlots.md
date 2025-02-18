# Generation of synteny plots of reference genomes and Hector's and Māui genome assemblies. 

## 1. Whole genome alingnments with minimap2
`Script for minimap2`
```
#!/bin/bash -e

#SBATCH --cpus-per-task 12
#SBATCH --job-name      alnnuc
#SBATCH --mem           50G
#SBATCH --time          06:00:00
#SBATCH --account       uoo02423
#SBATCH --output        %x_%j.out
#SBATCH --error         %x_%j.err
#SBATCH --hint          nomultithread

module purge
module load minimap2/2.28-GCC-12.3.0

minimap2 -t ${SLURM_CPUS_PER_TASK} \
-x asm10 tursiops_genome.fasta vaquita_genome.fasta > vaq_tur_aln.paf

minimap2 -t ${SLURM_CPUS_PER_TASK} \
-x asm5 maui_genome.fasta hectors_genome.fasta > maui_hec_aln.paf #asm5 for supspecies

minimap2 -t ${SLURM_CPUS_PER_TASK} \
-x asm10 vaquita_genome.fasta blue_genome.fasta > blue_vaquita_aln.paf
```
## 2. NGenomeSyn paf coordinates to link file
`Script for Paf2Link`
```
#!/bin/bash -e

#SBATCH --cpus-per-task 4
#SBATCH --job-name      pf2ln
#SBATCH --mem           10G
#SBATCH --time          01:00:00
#SBATCH --account       uoo02423
#SBATCH --output        %x_%j.out
#SBATCH --error         %x_%j.err
#SBATCH --hint          nomultithread

module purge
export PATH=/nesi/nobackup/uoo02423/Sebastian/SyntenyRev/NGenomeSyn/:$PATH

./bin/GetTwoGenomeSyn.pl Paf2Link vaq_tur_aln.paf 5000 vaq_tur_aln.link #5000 removes all alingments of fragments smaller than 5000 bp (default)

./bin/GetTwoGenomeSyn.pl Paf2Link maui_hec_aln.paf 5000 maui_hec_aln.link

./bin/GetTwoGenomeSyn.pl Paf2Link blue_vaquita_aln.paf 5000 blue_vaquita_aln.link
```
## 3. Generate synteny plots with NGenomeSyn. NGenomeSyn works with a .inf file. File example is mauihec.inf. In this file we specify the chromosome lenghts for each genome (.len files), and the .link file. Example file is "mauihec.inf"
```
################################### global parameters ############################################

SetParaFor = global

GenomeInfoFile1=mauichrgen.len
### Format (chr Start End ...other parameters)  Note: the order of chr present in the figure is the same with order defined in this file. if End Start was given in a line then reverse and complement of chr will be used
###  other parameters such as fill=red stroke-width=0  stroke=black stroke-opacity=1 fill-opacity=1 etc. Note: different lines could have different parameters


GenomeInfoFile2=hecchrgen.len  ##  GenomeInfoFile 2 means the 2th genome

LinkFileRef1VsRef2=maui_hec_aln.link
####### Format (chrA StartA EndA chrB StartB End ...other parameters)
### Note: link files could be occurred multiple times,as: LinkFileRef1VsRef2=


#Main = "main_Figure"  ##  the Fig Name  :MainRatioFontSize MainCor ShiftMainX  ShiftMainY 

############################### Figure ##########################################################


#############################  Canvas and image parameter configuration ##########################
#body=1200    ## default: 1200 and up/down/left/right) = (55,25,100,120); #CanvasHeightRitao=1.0 CanvasWidthRitao=1.0
#RotatePng=0  ## rotation angles for the figure





SetParaFor = Genome1  ## GenomeALL/GenomeX : setting parameter  for the ALL/1st genome

#ZoomChr=1.0          ## enlarge(>1) or reduce(<1) or equal (=1) the lenghth of chromosomes
#RotateChr=30         ## rotation angles for the chromosomes with clockwise(>0) or counterclockwise(<0)
#ShiftX=0
#ShiftY=0             ## move this genome
#MoveToX/Y

#ChrWidth=20          ## the height of chromosome
#LinkWidth=180        ## the link height between two genomes
#ChrSpacing=10        ## the space length between two chromosomes
#NormalizedScale=0    ## use the custome scale (=1) or the same with default genomes (=0)

ChrNameShow=1         ## show this genome chr name
ChrNameColor=green    ## Set Chr name color to green
#ChrNameShiftY=10      ## ChrName moves down 10
#ShowCoordinates=1     ## Show Coordinates . with other para [ScaleNum=10 ScaleUpDown ScaleUnit LabelUnit  LablefontsizeRatio  ]
##other rare use parameters, such as EndCurveRadian=3/ etc
#EndCurveRadian=1000
#ZoomChr=0.5
#RotateChr=30

SetParaFor = Genome2  ## GenomeALL/GenomeX : setting parameter for the ALL/2st genome
ChrNameShow=1
ChrNameColor=green
ChrNameShiftY=30
GenomeNameColor=blue
#ShowCoordinates=1     ## Show Coordinates . with other para [ScaleNum=10 ScaleUpDown ScaleUnit LabelUnit  LablefontsizeRatio  ]
SetParaFor =Link1     ## set parameter for the 1st occurrence of  LinkFileRef*VsRef* attributes
#StyleUpDown=line     ## change to straight line
```
`Script for NGenomeSyn`
```
#!/bin/bash -e

#SBATCH --cpus-per-task 4
#SBATCH --job-name      plosyn
#SBATCH --mem           10G
#SBATCH --time          01:00:00
#SBATCH --account       uoo02423
#SBATCH --output        %x_%j.out
#SBATCH --error         %x_%j.err
#SBATCH --hint          nomultithread

module purge
export PATH=/nesi/nobackup/uoo02423/Sebastian/SyntenyRev/NGenomeSyn/:$PATH

./bin/NGenomeSyn -InConf balmus_balric.inf -OutPut f_balmus_balric_synteny_mm2

./bin/NGenomeSyn -InConf bostau_balmus.inf -OutPut f_bostau_balmus_synteny_mm2

./bin/NGenomeSyn -InConf bostau_lagalb.inf -OutPut f_bostau_lagalb_synteny_mm2
```
