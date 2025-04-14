#!/bin/bash
#SBATCH --job-name=ArrayJob
#SBATCH --output=IQTREE_%A_%a.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01-00:00:00
#SBATCH --mem=30G
#SBATCH --qos=short

input=$1 # Input directory with alignments
output=$2 # Output directory for trees

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

# skip if done
if [[ -d $output/$name ]]; then
	echo "Done"
	exit 0
fi

mkdir -p $output/$name

snp-sites -v $INPUTFILE > $output/$name/$name.vcf

# Matriz de pairwise-dist
snp-dists -j 8 -m $INPUTFILE > $output/$name/$name.fasta.MD

#snps 
snp-sites $INPUTFILE -o $output/$name/$name.SNPs.fasta

#constant
fconstvar=$(snp-sites -C $INPUTFILE)

iqtree2 -s $output/$name/$name.SNPs.fasta \
	-m GTR+I+F+G4 \
	-nt 8 \
	--redo -fconst $fconstvar \
	--prefix $output/$name/$name

exit 0
