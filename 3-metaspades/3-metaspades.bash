#!/bin/bash

#SBATCH --job-name=metaspades
#SBATCH --account=dw30
#SBATCH --time=168:00:00
#SBATCH --partition=comp
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8192
#SBATCH --cpus-per-task=32
#SBATCH --qos=normal


### Time to run Metaspades! This is gonna take forever to finish running, so
### ideally set the walltime to the max and let it run overnight x 3. This
### script has the mem-per-cpu and cpus-per-task cranked way up to finish
### as fast as possible. Metaspades is used to assemble contigs and scaffolds.


# Module load
module load spades/3.13.1

# Vars
trimmomatic_outputs='/home/mnak0010/dw30/Michael/Projects/2021_may_metagenomeseq/2-trimmomatic/'
spades_path='/usr/local/spades/3.13.1/bin/spades.py'
output_dir='/home/mnak0010/dw30/Michael/Projects/2021_may_metagenomeseq/3-metaspades/outputs2/'
paired_suffix='paired.fq.gz'
unpaired_suffix='unpaired.fq.gz'
forward='1_paired.fq.gz'
reverse='2_paired.fq.gz'
threads='64'

# Main

# Sort out files in trimmomatic_outputs into paired and unpaired
mkdir -p ${trimmomatic_outputs}paired ${trimmomatic_outputs}unpaired
mv ${trimmomatic_outputs}*${unpaired_suffix} ${trimmomatic_outputs}unpaired/
mv ${trimmomatic_outputs}*${paired_suffix} ${trimmomatic_outputs}paired/

a=0
filelist=${trimmomatic_outputs}paired/*

for fl in $filelist
do
    if [ $a -eq 1 ]
    then
        a=0
    
    else
        a=1
        
        # Make the generic filename without the f or r designation
        mediumname1=${fl%_*}
	    mediumname2=${mediumname1%_*}_
        mediumname3=$(basename $mediumname2)
	    forwardname=${mediumname2}${forward}
        reversename=${mediumname2}${reverse}

	    mkdir ${output_dir}${mediumname3}

        $spades_path \
        --meta \
        -1 $forwardname \
        -2 $reversename \
        -t $threads \
        -o ${output_dir}${mediumname3}

    fi
done
