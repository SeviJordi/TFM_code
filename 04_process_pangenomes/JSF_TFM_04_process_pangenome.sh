#!/usr/bin/bash
#SBATCH --qos=short
#SBATCH --job-name=process_pangenome
#SBATCH -n 1
#SBATCH --mem=20GB
#SBATCH --output=%j_roary_proc.out 
#SBATCH --cpus-per-task=64 
#SBATCH --time=1-00:00:00

logthis () { #Log things
    echo $(date) "|" $1
}


# Jordi Sevilla Fortuny
# Script to process protein families alignments from core genome
# Usage: sbatch process_pangenome.sh [ST] [pan_output_folder] [Pangenome tool]

# Check if the user has provided the required arguments
if [ $# -ne 3 ]; then
    echo "Usage: sbatch process_pangenome.sh [ST] [pan_output_folder] [Pangenome tool]"
    exit 1
fi

# Set params depending on the tool used
if [[ "$3" == "panacota" ]]; then
    extension=".gen"
    aln_folder=$2/$(ls $2 | grep Align)
    files=($aln_folder/*.gen)

elif [[ "$3" == "panaroo" ]]; then
    extension=".aln.fas"
    aln_folder=$2/aligned_gene_sequences/
    files=($aln_folder/*.aln.fas)

elif [[ "$3" == "roary" ]]; then
    extension=".fa.aln"
    aln_folder=$2/pan_genome_sequences/
    gene_names=$(grep label $2/core_alignment_header.embl | cut -d"=" -f2)
    files=($(parallel echo $aln_folder/{}.fa.aln ::: $gene_names))

else
    echo "Pangenome tool not recognized. Please use 'panacota', 'panaroo' or 'roary'"
    exit 1
fi

# Functions to process the protein families alignments
mafft_align() { # $1: input file, $2: output directory, $3 extension
    
    name=$(basename $1 $3 | sed 's/current.//')
    mafft --thread 8 \
        --adjustdirection $1 |\
    sed -e 's/_[0-9]\{5\}\($\| .*\)//' \
        -e 's/;.*//g' \
        -e 's/_R_//g' > $2/$name.aln
}

export -f mafft_align

run_clipkit(){ # $1: input file, $2: output directory
    name=$(basename $1 .aln)
    clipkit $1 -o $2/$name.clipkit.aln
}

export -f run_clipkit

run_trimal(){ # $1: input file, $2: output directory
    name=$(basename $1 .aln)
    trimal -fasta -automated1 -in $1 -out $2/$name.trimal.aln
}

export -f run_trimal

### Main script ###
# Create output directories
mafft_dir=$2/mafft_alignments
clipkit_dir=$2/clipkit/aln
trimal_dir=$2/trimal/aln
treeshrink_dir=$2/treeshrink/
fixcore_dir=$2/fixcore

mkdir -p $mafft_dir $clipkit_dir $trimal_dir $treeshrink_dir $fixcore_dir

# Align protein families
logthis "Aligning protein families"
parallel -j 32 mafft_align {} $mafft_dir $extension ::: ${files[@]}

# Run clipkit
logthis "Running clipkit"
parallel -j 64 run_clipkit {} $clipkit_dir ::: $mafft_dir/*.aln

# Run trimal
logthis "Running trimal"
parallel -j 64 run_trimal {} $trimal_dir ::: $clipkit_dir/*.aln

# Run treeshrink
logthis "Preparing data for treeshrink"

mkdir -p $treeshrink_dir/data
for alignment in $mafft_dir/*.aln; do
    name=$(basename $alignment .aln)
    
    mkdir -p $treeshrink_dir/data/$name
    cp $alignment $treeshrink_dir/data/$name/alignment.fasta
    iqtree2 -s $alignment \
        -m GTR+F+I+G4 \
        -nt 2 \
        --prefix $treeshrink_dir/data/$name/tree

done

logthis "Executing treeshrink"
mkdir -p $treeshrink_dir/results

run_treeshrink.py -i $treeshrink_dir/data/ \
    -t tree.treefile -a alignment.fasta \
    -m "all-genes" \
    -o $treeshrink_dir/results > $treeshrink_dir/treeshrink.log


# Running hmmcleaner
logthis "Using hmmcleaner"

parallel -j 64 HmmCleaner.pl ::: $treeshrink_dir/data/**/alignment.fasta

# Fixcore
logthis "Fix core"
bash fixcore_scripts/fixcore_pipeline.sh $mafft_dir "${1}_${3}" $fixcore_dir

# generating whole lignments
logthis "Generating whole alignments"

mkdir $2/fixed_alignments

# Trimmal
AMAS.py concat -i $trimal_dir/* -f fasta -d dna -c 8 -o $2/fixed_alignments/${1}_${3}_trimal.aln

# Clipkit
AMAS.py concat -i $clipkit_dir/* -f fasta -d dna -c 8 -o $2/fixed_alignments/${1}_${3}_clipkit.aln

# treeshrink
AMAS.py concat -i $treeshrink_dir/results/**/output.fasta -f fasta -d dna -c 8 -o $2/fixed_alignments/${1}_${3}_treeshrink.aln

# Hmmcleaner
AMAS.py concat -i $treeshrink_dir/data/**/alignment_hmm.fasta -f fasta -d dna -c 8 -o $2/fixed_alignments/${1}_${3}_hmmcleaner.aln

# Fixcore
AMAS.py concat -i $fixcore_dir/DEF_GENES4/*.fasta -f fasta -d dna -c 8 -o $2/fixed_alignments/${1}_${3}_fixcore.aln

logthis "Done"
exit 0