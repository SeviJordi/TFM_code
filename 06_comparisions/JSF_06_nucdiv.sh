#!/bin/bash
#SBATCH --job-name=NucDiv
#SBATCH --output=nuc_div_%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --time=02-18:00:00
#SBATCH --mem=30G
#SBATCH --qos=medium

input=$1 # Input directory with alignments
output=$2 # Output directory for pairwise distances
prefix=$3 # Prefix for output files
# Check arguments
if [ "$#" -ne 3 ]; then
  echo "Incorrect number of arguments"
  echo "Usage: $0 [input] [output] [prefix]"
  exit 1
fi
# Check if input directory exists
if [ ! -d "$input" ]; then
  echo "Input directory does not exist: $input"
  exit 1
fi

# Create output
mkdir -p $output

# Analyze
mkdir -p $output/$prefix

calculate_div(){ # aln outdir
name=$(basename $1 .aln)
len=$(alnlen -i $1)

snp-dists -j 4 -m $1 | awk -v len=$len '{print $0"\t"len"\t" $3/len}' > $2/$name.DM 
}

export -f calculate_div

parallel -j 9 calculate_div ::: $input/* ::: $output/$prefix

exit 0
