#!/bin/bash
#SBATCH --job-name=nucdiv
#SBATCH --output=nucdiv_%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=01-00:00:00
#SBATCH --mem=25G
#SBATCH --qos=short

# Jordi Sevilla Fortuny
# Script to calculate nucleotide diversity from pairwise distances
# Usage: JSF_TFM_06_nucdiv.sh [input dir] [output file]
# Input dir: Directory with pairwise distance files
# Output file: File to save nucleotide diversity results
# Check arguments
if [ "$#" -ne 2 ]; then
  echo "Incorrect number of arguments"
  echo "Usage: $0 [input dir] [output file]"
  exit 1
fi
input_dir=$1 # Input directory with pairwise distance files
output_file=$2 # Output file for nucleotide diversity results
# Check if input directory exists
if [ ! -d "$input_dir" ]; then
  echo "Input directory does not exist: $input_dir"
  exit 1
fi

get_diversities(){
 	tool="raw"
        if [[ $1 =~ "snippy" ]]; then
                tool="snippy"
        fi
        
       	pan=$(basename $1 | cut -d"_" -f2)
        st=$(basename $1 | cut -d"_" -f1)
        for file in $1/*.DM; do 
                name=$(basename $file .DM)
                pi=$(awk -f nuc_div.awk $file)
                echo -e "$st\t$pan\t$tool\t$name\t$pi" >> $2 
        done
}

export -f get_diversities

parallel -j 64 get_diversities ::: $input_dir/* ::: $output_file
exit 0
