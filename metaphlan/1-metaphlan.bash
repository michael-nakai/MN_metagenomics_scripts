#!/bin/bash

#SBATCH --job-name=metaphlan
#SBATCH --account=dw30
#SBATCH -p comp
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64gb
#SBATCH --time=100:00:00

module load metaphlan/2.0

### Vars here
flopaired='/home/mnak0010/dw30_scratch/Michael/2021_may_metagenomeseq/2-trimmomatic/paired/flo/'
hamdipaired='/home/mnak0010/dw30_scratch/Michael/2021_may_metagenomeseq/2-trimmomatic/paired/hamdi/'

generalOutput='/home/mnak0010/dw30_scratch/Michael/2021_may_metagenomeseq/metaphlan/'
flooutFolder='/home/mnak0010/dw30_scratch/Michael/2021_may_metagenomeseq/metaphlan/flo/'
hamdioutFolder='/home/mnak0010/dw30_scratch/Michael/2021_may_metagenomeseq/metaphlan/hamdi/'

aggregateOutputs='/home/mnak0010/dw30_scratch/Michael/2021_may_metagenomeseq/metaphlan/aggregates/'
aggregatesFlo='/home/mnak0010/dw30_scratch/Michael/2021_may_metagenomeseq/metaphlan/aggregates/flo/'
aggregatesHamdi='/home/mnak0010/dw30_scratch/Michael/2021_may_metagenomeseq/metaphlan/aggregates/hamdi/'

metaphlan_database='/home/mnak0010/dw30_scratch/Michael/2021_may_metagenomeseq/scripts/humann2/metaphlan_database/'

cpu_number=12

### Main

mkdir -p $generalOutput $flooutFolder $hamdioutFolder $aggregateOutputs $aggregatesFlo $aggregatesHamdi

flo=$(ls "${flopaired}" | grep "1_paired.fq.gz")
hamdi=$(ls "${hamdipaired}" | grep "1_paired.fq.gz")

# Flo's stuff
for x in ${flo[@]}
do
        y=${x%_*}
        z=${y%_*}
        fw=${flopaired}${z}_1_paired.fq.gz
        rv=${flopaired}${z}_2_paired.fq.gz

        mkdir ${flooutFolder}${z}

        # Run metaphlan on the paired ends
        "/usr/local/metaphlan/2.0/MetaPhlAn/metaphlan2.py" \
            ${fw},${rv} \
            --bowtie2db $metaphlan_database \
            --bowtie2out ${flooutFolder}${z}'/'${z}'metagenome.bowtie2.bz2' \
            --nproc $cpu_number \
            --input_type fastq \
            -o ${flooutFolder}${z}'/'${z}'_profiled_metagenome.txt'
        
        # Copy the output to the aggregate folder
        cp ${flooutFolder}${z}'/'${z}'_profiled_metagenome.txt' ${aggregatesFlo}   
done

# Hamdi's stuff
for x in ${hamdi[@]}
do
        y=${x%_*}
        z=${y%_*}
        fw=${hamdipaired}${z}_1_paired.fq.gz
        rv=${hamdipaired}${z}_2_paired.fq.gz

        mkdir ${hamdioutFolder}${z}

        # Run metaphlan on the paired ends
        "/usr/local/metaphlan/2.0/MetaPhlAn/metaphlan2.py" \
            ${fw},${rv} \
            --bowtie2db $metaphlan_database \
            --bowtie2out ${hamdioutFolder}${z}'/'${z}'metagenome.bowtie2.bz2' \
            --nproc $cpu_number \
            --input_type fastq \
            -o ${hamdioutFolder}${z}'/'${z}'_profiled_metagenome.txt'
        
        # Copy the output to the aggregate folder
        cp ${hamdioutFolder}${z}'/'${z}'_profiled_metagenome.txt' ${aggregatesHamdi} 
done

# Merge results, then do heatmaps
"/usr/local/metaphlan/2.0/MetaPhlAn/utils/metaphlan_hclust_heatmap.py" \
    --in $flooutFolder'all_results_merged.txt' \
    --out $flooutFolder'all_results_heatmap.png' \
    -c bbcry \
    --top 25 \
    --minv 0.1 \
    -s log
"/usr/local/metaphlan/2.0/MetaPhlAn/utils/metaphlan_hclust_heatmap.py" \
    --in $hamdioutFolder'all_results_merged.txt' \
    --out $hamdioutFolder'all_results_heatmap.png' \
    -c bbcry \
    --top 25 \
    --minv 0.1 \
    -s log