#!/usr/bin/bash

#SBATCH --qos=short
#SBATCH --job-name=panaroo
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=40GB
#SBATCH --output=%j_panaroo.out
#SBATCH --cpus-per-task 32
#SBATCH --time=1-00:00:00 

# Script to run panaroo

# Usage: sbatch run_panaroo.sh annotations threshold data_dir

# Check number of arguments
if [ "$#" -ne 3 ]; then
    echo "Illegal number of parameters"
    echo "Usage: sbatch run_panaroo.sh annotations threshold data_dir"
    exit 1
fi

# Initialize values
annotations=$1
threshold=$2
data_dir=$3

# Move to data directory
cd $data_dir


# Execute panaroo in default mode
panaroo -i $annotations/*.gff \
    -o panaroo_default \
    --clean-mode sensitive \
    -t 32 \
    --alignment core \
    --aligner mafft \
    --core_threshold $threshold

# Execute panaroo in pretending to be panacota
panaroo -i $annotations/*.gff \
    -o panaroo_modifyed \
    --clean-mode sensitive \
    -t 32 \
    --alignment core \
    --aligner mafft \
    --core_threshold $threshold \
    -f 0.8 \
    --merge_paralogs


exit 0

