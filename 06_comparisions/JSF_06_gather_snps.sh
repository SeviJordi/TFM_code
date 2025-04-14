#!/usr/bin/bash

# Jordi Sevilla Fortuny
# Script to gather SNPs from different approaches
# Usage: JSF_TFM_06_gather_snps.sh [input dir] [output file]
# Input dir: Directory with VCF files
# Output file: File to save SNPs results

# Check arguments   
if [ "$#" -ne 2 ]; then
  echo "Incorrect number of arguments"
  echo "Usage: $0 [input dir] [output file]"
  exit 1
fi

input_dir=$1 # Input directory with VCF files
output_file=$2 # Output file for SNPs results

# Check if input directory exists
if [ ! -d "$input_dir" ]; then
  echo "Input directory does not exist: $input_dir"
  exit 1
fi

parse_file(){
    file=$1
    len=$(alnlen -i $file)
    python3 parse_vcf.py -i $file -l $len
}
export -f parse_file


find $input_dir -name \*.vcf -print0 | xargs -P 20 -n 1 -0 bash -c 'parse_file "$@"' {} >> $output_file

exit 0
