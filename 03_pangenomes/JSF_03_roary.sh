#!/usr/bin/bash

#SBATCH --qos=medium
#SBATCH --job-name=roary
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=40GB
#SBATCH --output=%j_roary 
#SBATCH --cpus-per-task 12
#SBATCH --time=1-00:00:00 

# Jordi Sevilla Fortuny
# Usage: sbatch roary.sh <annotations> <threshold> <data_dir>

# Check the number of arguments
if [ "$#" -ne 3 ]; then
    echo "Illegal number of parameters"
    echo "Usage: sbatch roary.sh <annotations> <threshold> <data_dir>"
    exit 1
fi

# Initialize the variables
annotations=$1
threshold=$2
data_dir=$3

# Move to data directory
cd $data_dir

# Execute roary in default mode
roary -e \
    --mafft \
    -p 12 \
    -i 80 \
    -f $data_dir/roary_default \
    -v -z \
    -cd $threshold $annotations/*.gff
