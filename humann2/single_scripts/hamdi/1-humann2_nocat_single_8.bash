#!/bin/bash

#SBATCH --job-name=HumanN2_hamLF3
#SBATCH --account=dw30
#SBATCH -p comp
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80gb
#SBATCH --time=100:00:00


### Load modules
module load bowtie2/2.3.5
module load humann/2.0
module load metaphlan/2.0
module load tbb/20180312oss

### Vars here
num=8

medium='/home/mnak0010/dw30_scratch/Michael/2021_may_metagenomeseq/2-trimmomatic/'

flopaired='/home/mnak0010/dw30_scratch/Michael/2021_may_metagenomeseq/2-trimmomatic/paired/flo/'
hamdipaired='/home/mnak0010/dw30_scratch/Michael/2021_may_metagenomeseq/2-trimmomatic/paired/hamdi/'

floCatFiles='/home/mnak0010/dw30_scratch/Michael/2021_may_metagenomeseq/humann2/catfiles/flo/'
hamdiCatFiles='/home/mnak0010/dw30_scratch/Michael/2021_may_metagenomeseq/humann2/catfiles/hamdi/'

flooutFolder='/home/mnak0010/dw30_scratch/Michael/2021_may_metagenomeseq/humann2/flo2/'
hamdioutFolder='/home/mnak0010/dw30_scratch/Michael/2021_may_metagenomeseq/humann2/hamdi2/'

uniref_database='/home/mnak0010/dw30_scratch/Michael/shared_databases/uniref'
chocophlan_database='/home/mnak0010/dw30_scratch/Michael/shared_databases/chocophlan'

metaphlan_database='/home/mnak0010/dw30_scratch/Michael/2021_may_metagenomeseq/scripts/humann2/metaphlan_database/'


### MAIN
mkdir -p $flooutFolder $hamdioutFolder $floCatFiles $hamdiCatFiles

hamdi=(TG-HF1-LDE4414_L3_1_paired.fq.gz TG-HF2-LDE4415_L3_1_paired.fq.gz TG-HF3-LDE4416_L3_1_paired.fq.gz TG-HF4-LDE4417_L3_1_paired.fq.gz TG-HF5-LDE4418_L3_1_paired.fq.gz \
        TG-HF6-LDE4419_L3_1_paired.fq.gz TG-LF1-LDE4420_L3_1_paired.fq.gz TG-LF2-LDE4421_L3_1_paired.fq.gz TG-LF3-LDE4422_L3_1_paired.fq.gz TG-LF4-LDE4423_L3_1_paired.fq.gz \
        TG-LF5-LDE4424_L3_1_paired.fq.gz TG-LF6-LDE4425_L3_1_paired.fq.gz)

x=${hamdi[$num]}
y=${x%_*}
z=${y%_*}
fw=${hamdipaired}${z}_1_paired.fq.gz
rv=${hamdipaired}${z}_2_paired.fq.gz
catname=${z}"_cat.fq.gz"

# Cat the paired end files together (https://github.com/biobakery/humann#humann-30-and-paired-end-sequencing-data)
#cat $fw $rv > ${floCatFiles}${catname}

#Run HumanN2 on the catfile
"/usr/local/humann/2.0/humannv2.0/bin/humann2" \
        --input ${hamdiCatFiles}${catname} \
        --output ${hamdioutFolder}${z} \
        --protein-database $uniref_database \
        --nucleotide-database $chocophlan_database \
        --metaphlan-options "--bowtie2db $metaphlan_database"