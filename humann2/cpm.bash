#!/bin/bash

#SBATCH --job-name=cpm
#SBATCH --account=dw30
#SBATCH -p comp
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80gb
#SBATCH --time=100:00:00


### Load modules
module load humann/2.0
module load bowtie2/2.3.5
module load tbb/20180312oss

### Vars here
pathcov=/home/mnak0010/dw30_scratch/Michael/2021_may_metagenomeseq/humann2/all_tsvs/pathabun/*

pathcout='/home/mnak0010/dw30_scratch/Michael/2021_may_metagenomeseq/humann2/outputs/all_tsvs/pathabun/'


### MAIN
mkdir -p $geneout $pathcout $pathaout

for x in ${pathcov[@]}
do
        c=${x%.*}
        filename=${c##*/}
        
        # Change RPK to CPM
        "/usr/local/humann/2.0/humannv2.0/bin/humann2_renorm_table" \
                --input $x \
                --output ${pathcout}${filename}'_cpm.tsv' \
                --units cpm \
                --update-snames
done
