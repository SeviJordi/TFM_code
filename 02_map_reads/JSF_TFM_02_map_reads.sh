#!/usr/bin/bash

#SBATCH --partition=desa
#SBATCH --job-name=snippy
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=16030
#SBATCH --output=%j_download_sra.out
#SBATCH --cpus-per-task 30

# Jordi Sevilla Fortuny
# Script to map reads of a certain ST to a refference
# Usage: JSF_TFM_02_map_reads.sh [reads_dir] [ref] [path to data]

# Check args:
if [ "$#" -ne 3 ]; then
  echo "Incorrect number of arguments"
  exit 1
fi

# Initialize variables
reads_dir=$1
ref=$2
outdir=$3

# create snippy_dir
mkdir -p $outdir/snippy

#Function to map
map_snippy(){

    snippy --outdir $3/$1 \
        --ref $4 \
        --R1 $2/$1/${1}_1.fq \
        --R2 $2/$1/${1}_2.fq
}

export -f map_snippy

# process
parallel -j 10 map_snippy ::: $(ls $reads_dir) ::: $reads_dir ::: $snippy_dir ::: $ref

# Gather all results with snippy core

real_snippy_path=$(realpath $outdir/snippy)
real_ref_path=$(realpath $ref)

mkdir -p $outdir/snippy_core
cd $outdir/snippy_core

snippy_core --ref $real_ref_path $real_snippy_path/*

exit 
