#!/bin/bash

#SBATCH --job-name=HumanN2_hamdi1
#SBATCH --account=dw30
#SBATCH -p comp
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80gb
#SBATCH --time=168:00:00


### Load modules
module load humann/2.0
module load bowtie2/2.3.5
module load tbb/20180312oss

### Vars here
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

flo=$(ls "${flopaired}" | grep "1_paired.fq.gz")
hamdi=(TG-HF3-LDE4416_L3_1_paired.fq.gz TG-HF4-LDE4417_L3_1_paired.fq.gz TG-HF5-LDE4418_L3_1_paired.fq.gz)

# Hamdi's stuff
for x in ${hamdi[@]}
do
        y=${x%_*}
        z=${y%_*}
        fw=${hamdipaired}${z}_1_paired.fq.gz
        rv=${hamdipaired}${z}_2_paired.fq.gz
        catname=${z}"_cat.fq.gz"

        # Cat the paired end files together (https://github.com/biobakery/humann#humann-30-and-paired-end-sequencing-data)
        #cat $fw $rv > ${hamdiCatFiles}${catname}

        #Run HumanN2 on the catfile
        "/usr/local/humann/2.0/humannv2.0/bin/humann2" \
                --input ${hamdiCatFiles}${catname} \
                --output ${hamdioutFolder}${z} \
                --protein-database $uniref_database \
                --nucleotide-database $chocophlan_database \
                --metaphlan-options "--bowtie2db $metaphlan_database"
done
