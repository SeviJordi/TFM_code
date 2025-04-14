#!/bin/bash
#SBATCH --job-name=ArrayJob
#SBATCH --output=AMAS_%A_%a.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=60G
#SBATCH --qos=short

input=$1 # Input directory with alignments
output=$2 # Output directory for summaries

# Check arguments
if [ "$#" -ne 2 ]; then
  echo "Incorrect number of arguments"
  echo "Usage: $0 [input] [output]"
  exit 1
fi

# Create output
mkdir -p $output

# Set array
FILES=($input/*)
INPUTFILE=${FILES[$SLURM_ARRAY_TASK_ID]}

# Analyze
name=$(basename $INPUTFILE .aln)

if [[ ! -f  $output/$name.summary ]]; then
	AMAS.py summary -i $INPUTFILE -f fasta -d dna -o $output/$name.summary
fi

exit 0
