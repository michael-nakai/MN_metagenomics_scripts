#!/bin/bash

#SBATCH --job-name=trimmomatic
#SBATCH --account=dw30
#SBATCH --time=168:00:00
#SBATCH --partition=comp
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8192
#SBATCH --cpus-per-task=16
#SBATCH --qos=normal


### Overview
### Trimmomatic is used to comb through the fastq files and:
###     -Trim Illumina adaptors
###     -Remove leading/trailing low quality bases
###     -Pair forward and reverse reads
###     -Cut

# Module loading
module load trimmomatic/0.38

# Vars
trimmomatic_path='/usr/local/trimmomatic/0.38/trimmomatic-0.38.jar'
filelist=/home/mnak0010/dw30/Michael/Projects/2021_may_metagenomeseq/Metagenome_Seq_May2021/raw-files/*
forward='1.fq.gz'
reverse='2.fq.gz'
outputfolder='/home/mnak0010/dw30/Michael/Projects/2021_may_metagenomeseq/2-trimmomatic/'

# Trimmomatic settings
illumina_adaptors='TruSeq3-PE.fa:2:30:10'
leading='3' #Remove leading low quality or Ns below quality x
trailing='3' #Remove trailing low quality or Ns below quality x
sliding_window='4:15' #x-base wide sliding window, cutting when the average quality per base drops below y
minimum_read_length=36

# Main
a=0

for fl in $filelist
do
    if [ $a -eq 1 ]
    then
        a=0
    
    else
        a=1
        
        # Make the generic filename without the f or r designation
        mediumname=${fl%_*}_
        forwardname=${mediumname}${forward}
        reversename=${mediumname}${reverse}
        o_f_paired=${outputfolder}${mediumname##*/}${forward%%.*}_paired.fq.gz
        o_f_unpaired=${outputfolder}${mediumname##*/}${forward%%.*}_unpaired.fq.gz
        o_r_paired=${outputfolder}${mediumname##*/}${reverse%%.*}_paired.fq.gz
        o_r_unpaired=${outputfolder}${mediumname##*/}${reverse%%.*}_unpaired.fq.gz

        java -jar $trimmomatic_path PE \
	-phred33 \
        $forwardname \
        $reversename \
        $o_f_paired \
        $o_f_unpaired \
        $o_r_paired \
        $o_r_unpaired \
        ILLUMINACLIP:$illumina_adaptors \
        LEADING:$leading \
        TRAILING:$trailing \
        SLIDINGWINDOW:$sliding_window \
        MINLEN:$minimum_read_length

    fi
done


