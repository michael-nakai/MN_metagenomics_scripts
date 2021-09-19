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
genefam=/fs03/dw30/Michael/2021_may_metagenomeseq/humann2/outputs/genefam-cpm/*

geneout='/home/mnak0010/dw30_scratch/Michael/2021_may_metagenomeseq/humann2/outputs/all_tsvs/genefam_EC_regrouped/'


### MAIN
mkdir -p $geneout

for x in ${genefam[@]}
do
        c=${x%.*}
        filename=${c##*/}
        
        # Regroup according to EC
        "/usr/local/humann/2.0/humannv2.0/bin/humann2_regroup_table" \
                --input $x \
                --output ${pathcout}${filename}'_cpm_regrouped.tsv' \
                --custom '/home/mnak0010/dw30_scratch/Michael/2021_may_metagenomeseq/humann2/utility_db/utility_mapping/map_level4ec_uniref90.txt.gz'
done
