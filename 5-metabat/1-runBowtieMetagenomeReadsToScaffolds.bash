#!/bin/bash

#SBATCH --job-name=ReadsToScaffolds
#SBATCH --account=dw30
#SBATCH -p comp
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64gb
#SBATCH --time=100:00:00


### Vars here
medium='/home/mnak0010/dw30/Michael/Projects/2021_may_metagenomeseq/2-trimmomatic/'

flopaired='/home/mnak0010/dw30/Michael/Projects/2021_may_metagenomeseq/2-trimmomatic/paired/flo/'
flounpaired='/home/mnak0010/dw30/Michael/Projects/2021_may_metagenomeseq/2-trimmomatic/unpaired/flo/'

hamdipaired='/home/mnak0010/dw30/Michael/Projects/2021_may_metagenomeseq/2-trimmomatic/paired/hamdi/'
hamdiunpaired='/home/mnak0010/dw30/Michael/Projects/2021_may_metagenomeseq/2-trimmomatic/unpaired/hamdi/'

scaffoldFolder='/home/mnak0010/dw30/Michael/Projects/2021_may_metagenomeseq/3-metaspades/outputs/'

flooutFolder='/home/mnak0010/dw30/Michael/Projects/2021_may_metagenomeseq/4-metabat/outputs/flo/'
hamdioutFolder='/home/mnak0010/dw30/Michael/Projects/2021_may_metagenomeseq/4-metabat/outputs/hamdi/'


### MAIN
mkdir -p $flooutFolder $hamdioutFolder

flo=$(ls "${medium}paired/flo/" | grep "1_paired.fq.gz")
hamdi=$(ls "${medium}paired/hamdi/" | grep "1_paired.fq.gz")

module load bowtie2/2.3.5
module load samtools/1.9
module load bedtools/2.26.0

# Flo's stuff
for x in ${flo[@]}
do
        y=${x%_*}
        z=${y%_*}
        i=${z}_1_paired.fq.gz
        fw=${flopaired}${z}_1_paired.fq.gz
        rv=${flopaired}${z}_2_paired.fq.gz
        u1=${flounpaired}${z}_1_unpaired.fq.gz
        u2=${flounpaired}${z}_2_unpaired.fq.gz
        scaffold=${scaffoldFolder}${z}'_/scaffolds.fasta'
        bowtie2 -p 12 -x $scaffold -1 $fw -2 $rv -U $u1,$u2 -S ${flooutFolder}${z}.sam >& ${flooutFolder}${z}.log
        samtools view -S -b ${flooutFolder}${z}.sam > ${flooutFolder}${z}.bam
        samtools sort ${flooutFolder}${z}.bam -o ${flooutFolder}${z}.sorted.bam
        rm -f ${flooutFolder}${z}.sam
done

# Hamdi's stuff
for x in ${hamdi[@]}
do
        y=${x%_*}
        z=${y%_*}
        i=${z}_1_paired.fq.gz
        fw=${hamdipaired}${z}_1_paired.fq.gz
        rv=${hamdipaired}${z}_2_paired.fq.gz
        u1=${hamdiunpaired}${z}_1_unpaired.fq.gz
        u2=${hamdiunpaired}${z}_2_unpaired.fq.gz
        scaffold=${scaffoldFolder}${z}'_/scaffolds.fasta'
        bowtie2 -p 12 -x $scaffold -1 $fw -2 $rv -U $u1,$u2 -S ${hamdioutFolder}${z}.sam >& ${hamdioutFolder}${z}.log
        samtools view -S -b ${hamdioutFolder}${z}.sam > ${hamdioutFolder}${z}.bam
        samtools sort ${hamdioutFolder}${z}.bam -o ${hamdioutFolder}${z}.sorted.bam
        rm -f ${flooutFolder}${z}.sam
done
