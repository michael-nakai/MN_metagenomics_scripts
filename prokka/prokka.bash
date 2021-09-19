#!/bin/bash

#SBATCH --job-name=prokka_scaffolds
#SBATCH --account=dw30
#SBATCH -p comp
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80gb
#SBATCH --time=168:00:00

### Load modules
module load prokka/1.14.6

### Vars here
output_dir='/home/mnak0010/dw30_scratch/Michael/2021_may_metagenomeseq/prokka/'

### Main
mkdir -p output_dir

for scaffold in ~/dw30_scratch/Michael/2021_may_metagenomeseq/3-metaspades/outputs/*/scaffolds.fasta
do
    mid=${scaffold%-*}
    samplename=${mid##*/}
    
    prokka \
        --force \
        --outdir ${output_dir}${samplename} \
        --prefix $samplename \
        --centre X --compliant \
	$scaffold
done
