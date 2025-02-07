#!/usr/bin/bash

#SBATCH --partition=prod
#SBATCH --job-name=download_reads
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=16030
#SBATCH --output=%j_download_sra.out
#SBATCH --cpus-per-task 96

# Jordi Sevilla Fortuny
# Script to download SRA reads
# The script takes as input the output of the JSF_TFM_01_biosample2sra.sh script 
#  for one ST, the ST and the path where the reads will be stored.

# Usage: JSF_TFM_01_download_reads.sh [accecions file] [ST] [path to data]

# Check args
if [ "$#" -ne 3 ]; then
  echo "Incorrect number of arguments"
  exit 1
fi

# Initialize 
data=$3
st=$2
accesions=$1

# Go to the working directory and create reads folder
cd $data
mkdir -p $st/reads
cd $st/reads

# Function to download
download_sra() { # [args]

echo $1 > .temp_$2.txt

while read name biosample SRAs; do

    if [[ -n $SRAs ]]; then
        mkdir $name

        fasterq-dump -O $name -e 8 --split-files $SRAs

        # concatenate reads if there are more than one experiment
        cat $name/*1.fastq > $name/${name}_1.fq
        cat $name/*2.fastq > $name/${name}_2.fq

        echo "Done for $name"
    fi
done < .temp_$2.txt

rm .temp_$2.txt

}

export -f download_sra

# Download
parallel -j 12 -a $accesions  download_sra {} {#}

exit 0