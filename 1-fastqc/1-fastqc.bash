#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --account=dw30
#SBATCH --time=168:00:00
#SBATCH --partition=comp
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4096
#SBATCH --cpus-per-task=8
#SBATCH --qos=normal


### Overview
### Just standard fastqc


module load fastqc/0.11.9

rawdata=/home/mnak0010/dw30/Michael/Projects/2021_may_metagenomeseq/Metagenome_Seq_May2021/raw-files/*
outputdir=/home/mnak0010/dw30/Michael/Projects/2021_may_metagenomeseq/1-fastqc/

fastqc ${rawdata[@]} --outdir $outputdir
