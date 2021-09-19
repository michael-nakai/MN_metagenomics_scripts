#!/bin/bash

#SBATCH --job-name=busco
#SBATCH --account=dw30
#SBATCH -p comp
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH --time=100:00:00

bin_folders='~/dw30_scratch/Michael/2021_may_metagenomeseq/4-metabat/metabat_outputs/outputs1/' # Unused for now, not sure why for loop doesnt work with this var
database='/home/mnak0010/dw30_scratch/Michael/shared_databases/busco/'
general_output='/home/mnak0010/dw30_scratch/Michael/2021_may_metagenomeseq/5-busco/outputs/'


### Main

source ~/dw30_scratch/Michael/miniconda/main/miniconda/bin/activate
conda activate /home/mnak0010/dw30_scratch/Michael/conda_envs/busco

mkdir -p $general_output
cd $general_output

lineages=('bacteria_odb10' 'eukaryota_odb10' 'archaea_odb10')

for fld in ~/dw30_scratch/Michael/2021_may_metagenomeseq/4-metabat/metabat_outputs/outputs1/*/*/
do
    mid=${fld%*/*/}
    samplename_full=${mid##*/}
    samplename=${samplename_full%-*}

    for lineage in ${lineages[@]}
    do
        busco -i $fld \
            -o ${samplename}_${lineage} \
            -m genome \
            --lineage $lineage \
            --offline \
            --download_path $database
    done

done
