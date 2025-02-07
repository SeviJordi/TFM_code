#!/usr/bin/bash
# Jordi Sevilla Fortuny

# Script to find sra ID from biosample ID

# The input for this script is a tab delimited file 
#  with the sample names as the first column and the biosample id as the second column.

# The output is a tab separated file with the sample name and the biosample id as the first and 
#  second columns and then, all the sra IDs found separated by tabs.

# usage:
# bash JSF_TFM_01_biosample2sra.sh [input file] [output file]

# Check arguments
if [ "$#" -ne 2 ]; then
  echo "Incorrect number of arguments"
  exit 1
fi

# Finction to search 
get_sra(){

    query=$(echo "$2[BioSample] AND \
        "library layout paired"[Properties] AND \
        "platform illumina"[Properties]")
    echo $query
    acc=$(esearch -db sra -query "$query" </dev/null |\
            efetch -format docsum |\
            xtract -pattern Runs -element Run@acc -sep " " -tab " " |
            awk '{printf $1" "}')
    
    echo "$1 $2 $acc" >> $3
}

export -f get_sra

while read name sample; do
    get_sra $name $sample $2
done < $1


exit 0