#!/usr/bin/bash

#SBATCH --qos=short
#SBATCH --job-name=snippy_subset
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=40GB
#SBATCH --output=%j_subset_snippy.out
#SBATCH --cpus-per-task 32 
#SBATCH --time=1-00:00:00

# Jordi Sevilla Fortuny
# Script pto subset core alignments with the intersection with the mapping
# Uso: JSF_TFM_03_subset_snippy.sh [referencia] [path to alignments] [snippy alignment] [st] [prefix]

if [ "$#" -ne 5 ]; then
  echo "Incorrect number of arguments"
  echo "Usage: JSF_TFM_05_subset_pangenomes.sh [reference] [path to alignments] [snippy alignment] [st] [prefix]"
  exit 1
fi

ref=$1
snippy_aln=$3
aln_dir=$2
st=$4
prefix=$5

files_dir=$(dirname $snippy_aln)
outdir=$files_dir/$prefix
mkdir -p $outdir

# Create multifasta with cnsensus
rm $outdir/cons_genes.fasta

do_consensus(){

	name=$(basename $1 .aln)
	consensus.py -n $name < $1 >> $2/cons_genes.fasta
}

export -f do_consensus

parallel -k -j 32 do_consensus ::: $aln_dir/* ::: $outdir


# Remove gaps from consensus
sed -i 's/-//g' $outdir/cons_genes.fasta

# Blast against ref
blastn -query $outdir/cons_genes.fasta \
       -subject $ref \
       -out $outdir/blast_results.tsv \
       -outfmt 6 \
       -perc_identity 80 \
       -qcov_hsp_perc 99

# filter blast results
awk -v column=1 -f remove_dups.awk $outdir/blast_results.tsv > $outdir/blast_results_filtered.tsv

# Generate list of positions to subset
regions=$(awk 'BEGIN{ORS=", "}{if ($9 < $10) {printf $9 "-" $10 ORS;}\
                            if ($9 > $10) {printf $10 "-" $9 ORS;}}'\
                            $outdir/blast_results.tsv | sed 's/..$//')


# subset alignment
extractalign -sequence $snippy_aln \
    -regions "$regions" \
    -outseq $outdir/${prefix}_snippy.aln

# store gene families from mapping
## Create partitions file
awk -f blast2partitions.awk $outdir/blast_results_filtered.tsv > $outdir/${prefix}_snippy.partitions

mkdir $outdir/snippy_genes
curdir=$(pwd)
realref=$(realpath $ref)

cd $outdir/snippy_genes
AMAS.py split -i $realref -d dna -f fasta -l ../${prefix}_snippy.partitions -j

for file in *out.fas; do
	newname=$(echo $file | sed  -e 's/^_//g' -e 's/-out.fas/.aln/g')
	mv $file $newname
done
cd $curdir

# Create a gene list
cut -f1 $outdir/blast_results_filtered.tsv > $outdir/$prefix.genes

# For panacota, genes have to be changued
if [[ $prefix =~ "panacota" ]]; then
	sed -i 's/'$st'/'$st'-/g' $outdir/$prefix.genes
fi

mkdir -p $outdir/raw_genes
rm $outdir/raw_genes/*
 
for gene in $(cat $outdir/$prefix.genes); do
	cp $aln_dir/$gene.aln $outdir/raw_genes
done

AMAS.py concat -i $outdir/raw_genes/* -f fasta -d dna -c 8 -t $outdir/${prefix}_subset.aln -p $outdir/${prefix}_subset.partitions
 

# Subset corrections
path_to_cor=$(realpath $aln_dir/..)


## clipkit
mkdir -p $outdir/clipkit_genes
rm $outdir/clipkit_genes/*

for gene in $(cat $outdir/$prefix.genes); do

        if [[ -f $path_to_cor/clipkit/aln/$gene.clipkit.aln ]]; then
		cp $path_to_cor/clipkit/aln/$gene.clipkit.aln $outdir/clipkit_genes/$gene.aln
	fi
done

AMAS.py concat -i $outdir/clipkit_genes/* -f fasta -d dna -c 8 -t $outdir/${prefix}_clipkit.aln -p $outdir/${prefix}_clipkit.partitions

# TRIMAI
mkdir -p $outdir/trimal_genes
rm $outdir/trimal_genes/*
for gene in $(cat $outdir/$prefix.genes); do

        cp $path_to_cor/trimal/aln/$gene.trimal.aln $outdir/trimal_genes/$gene.aln
done

AMAS.py concat -i $outdir/trimal_genes/* -f fasta -d dna -c 8 -t $outdir/${prefix}_trimal.aln -p $outdir/${prefix}_trimal.partitions

# Treeshrink
mkdir -p $outdir/treeshrink_genes
rm $outdir/treeshrink_genes/*

for gene in $(cat $outdir/$prefix.genes); do

        cp $path_to_cor/treeshrink/results/$gene/output.fasta $outdir/treeshrink_genes/$gene.aln
done


AMAS.py concat -i $outdir/treeshrink_genes/* -f fasta -d dna -c 8 -t $outdir/${prefix}_treeshrink.aln -p $outdir/${prefix}_treeshrink.partitions

# hmmcleaner
mkdir -p $outdir/hmmcleaner_genes
rm $outdir/hmmcleaner_genes/*

for gene in $(cat $outdir/$prefix.genes); do

        cp $path_to_cor/treeshrink/data/$gene/alignment_hmm.fasta $outdir/hmmcleaner_genes/$gene.aln
done

AMAS.py concat -i $outdir/hmmcleaner_genes/* -f fasta -d dna -c 8 -t $outdir/${prefix}_hmmcleaner.aln -p $outdir/${prefix}_hmmcleaner.partitions


# fixcore
mkdir -p $outdir/fixcore_genes
rm $outdir/fixcore_genes/*

for gene in $(cat $outdir/$prefix.genes); do
	cp $path_to_cor/fixcore/DEF_GENES4/$gene.*fasta $outdir/fixcore_genes/$gene.aln

done

AMAS.py concat -i $outdir/fixcore_genes/* -f fasta -d dna -c 8 -t $outdir/${prefix}_fixcore.aln -p $outdir/${prefix}_fixcore.partitions

# end


