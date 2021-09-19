#!/bin/bash

#SBATCH --job-name=bowtie2build
#SBATCH --account=dw30
#SBATCH -p comp
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64gb
#SBATCH --time=100:00:00


### Overview
### This script makes a bowtie2 index for each sample cohort you have.
### For example, this cohort had Hamdi and Flo's samples, so I made a 
### separate index for Hamdi and Flo using one of their samples.

### VARS
module load bowtie2/2.3.5

medium='/home/mnak0010/dw30_scratch/Michael/2021_may_metagenomeseq/2-trimmomatic/'

outfolder='/home/mnak0010/dw30_scratch/Michael/2021_may_metagenomeseq/scripts/bowtie2-build'

indeces='/home/mnak0010/dw30_scratch/Michael/2021_may_metagenomeseq/3-metaspades/outputs/'


### MAIN
cd $outfolder

flo=$(ls "${medium}paired/flo/" | grep "1_paired.fq.gz")
hamdi=$(ls "${medium}paired/hamdi/" | grep "1_paired.fq.gz")

for x in ${flo[@]}
do
        y=${x%_*}
        z=${y%_*}
        indecesFull=${indeces}${z}'_/scaffolds.fasta'
        bowtie2-build ${indecesFull} ${indecesFull}
done

for x in ${hamdi[@]}
do
        y=${x%_*}
        z=${y%_*}
        indecesFull=${indeces}${z}'_/scaffolds.fasta'
        bowtie2-build ${indecesFull} ${indecesFull}
done

