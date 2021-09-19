#!/bin/bash

#SBATCH --job-name=metaphlan_heatmap
#SBATCH --account=dw30
#SBATCH -p comp
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
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