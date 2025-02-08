#!/usr/bin/bash

#SBATCH --partition=prod
#SBATCH --job-name=panacota
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=16030
#SBATCH --output=%j_panacota
#SBATCH --cpus-per-task 32

# Jordi Sevilla Fortuny
# Script to run panacota
# Usage: ./run_panacota.sh ST genomes threshold output db_genomes
# The genomes file is a file with one filename per line.
# The db_genomes folder contains the genomes in fasta format.

# Check number of arguments
if [ "$#" -ne 5 ]; then
    echo "Illegal number of parameters"
    echo "Usage: ./run_panacota.sh ST genomes threshold output db_genomes"
    exit 1
fi

# Function to log things
logthis () { #Log things
    echo $(date) "|" $1
}


# Inicializar
ST=$1
genomes=$2
threshold=$3
output=$4
output_dir=$output/panacota
db_genomes=$5
mkdir -p $output_dir


# Annotate the genomes
mkdir -p $output_dir/annotation

logthis "starting annotation step"

run_panacota.py annotate --nbcont 2000 \
    --l90 500 \
    -n KLPN \
    --cutn 0 \
    -r $output_dir/annotation \
    --threads 32 \
    -l $genomes \
    -q \
    -d $db_genomes

# Starting pangenome step

logthis "starting pangenome step"
mkdir -p $output_dir/pangenome

run_panacota.py pangenome -l $output_dir/annotation/LSTINFO-*.lst \
    -n "KLPN_$ST" \
    -d $output_dir/annotation/Proteins \
    -o $output_dir/pangenome \
    --threads 32 \
    -q 


# Starting corepears step

logthis "starting corepers step"

mkdir -p $output_dir/corepers

run_panacota.py corepers -p $output_dir/pangenome/Pan*lst \
    -t $threshold \
    -o $output_dir/corepers


# Starting core step

logthis "starting core step"
mkdir -p $output_dir/alignment

run_panacota.py align --threads 32 \
    -F \
    -n "KLPN_$ST" \
    -c $output_dir/corepers/PersGenome*$threshold*.lst \
    -l $output_dir/annotation/LSTINFO-*.lst \
    -d $output_dir/annotation \
    -o $output_dir/alignment

# Fix sample names
for file in $output_dir/alignment/Align-KLPN_$ST/*.gen; do 
    sed -i 's/>\([^\.]*\)\.\([^\.]*\)\.\([^\.]*\)\..*/>\1_\2_\3/g' $file
 done

logthis "Panacota done"
