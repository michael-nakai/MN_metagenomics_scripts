#!/bin/bash

#SBATCH --job-name=checkm
#SBATCH --account=dw30
#SBATCH -p comp
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64gb
#SBATCH --time=100:00:00

module load checkm/1.1.3

### VARS
medium='/home/mnak0010/dw30_scratch/Michael/2021_may_metagenomeseq/4-metabat/metabat_outputs/outputs2/'

outFolder='/home/mnak0010/dw30_scratch/Michael/2021_may_metagenomeseq/4-metabat/checkm_outputs/outputs2/'

numOfThreads=16


### MAIN
mkdir -p $outFolder

folderNames=(${medium}*)
fls=($(ls "${medium}"* | grep "metabat-bins16"))
i=0

for binFolder in ${fls[@]}
do
    echo $i
    echo '---------------------------------------------'
    finalOutFolder=${outFolder}$(basename ${folderNames[$i]})
    mkdir $finalOutFolder
    checkm lineage_wf --tab_table -x fa -t $numOfThreads ${folderNames[$i]}'/'${binFolder} $finalOutFolder
    let "i+=1"
done