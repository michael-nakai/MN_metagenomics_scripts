#!/bin/bash

#SBATCH --job-name=metabat
#SBATCH --account=dw30
#SBATCH -p comp
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64gb
#SBATCH --time=167:00:00

module load metabat/2.15.5

### VARS
medium='/home/mnak0010/dw30_scratch/Michael/2021_may_metagenomeseq/4-metabat/outputs/'
floBamFolder="${medium}flo/"
hamdiBamFolder="${medium}hamdi2/"

scaffoldFolder='/home/mnak0010/dw30_scratch/Michael/2021_may_metagenomeseq/3-metaspades/outputs2/'

mainoutFolder='/home/mnak0010/dw30_scratch/Michael/2021_may_metagenomeseq/4-metabat/metabat_outputs/'
hamdioutFolder='/home/mnak0010/dw30_scratch/Michael/2021_may_metagenomeseq/4-metabat/metabat_outputs/outputs2/'

### MAIN
mkdir -p $mainoutFolder $hamdioutFolder

hamdi=("TG-LF2-LDE4421_L3.sorted.bam" "TG-LF6-LDE4425_L3.sorted.bam")

for sortedBam in ${hamdi[@]}
do
    cd $hamdioutFolder
    y=${sortedBam%_*}
    genname=${y}'_L3' # IMPORTANT: These samples ended with _L3, but this should be changed per cohort if needed
    mkdir $genname
    cd $genname
    contig=${scaffoldFolder}${genname}'_/scaffolds.fasta' # Same for here, check the leading _
    fullsortedBam=${hamdiBamFolder}${sortedBam}
    runMetaBat.sh -m 2000 -t 16 $contig $fullsortedBam
done