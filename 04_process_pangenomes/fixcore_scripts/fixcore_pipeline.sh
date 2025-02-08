#!/bin/bash

# Adaptation of Neris García González pipeline to fix core alignments
# Usage: ./fix_core_alignments.sh <alignments> <prefix> <pathout>

logthis () { #Log things
    echo $(date) "|" $1
}

# Check arguments
if [ $# -ne 3 ]; then
    echo "Usage: ./fix_core_alignments.sh <alignments> <prefix> <pathout>"
    exit 1
fi

# Init variables
ALIGNMENTS=$(echo $1)
PREFIX=$(echo $2)
PATHOUT=$(echo $3)
dir=$(dirname $0)

mkdir -p $PATHOUT

# Set working tree
PATHFIXCORE=$PATHOUT/FIXCORE4
PATHCONS="$PATHFIXCORE/CONS"
PATHALG="$PATHFIXCORE/ALN"
PATHGENETREE="$PATHFIXCORE/GENETREES"
PATHAMAS="$PATHOUT/DEF_GENES4"
PATHTREE="$PATHOUT/Phylo4-$PREFIX"
TMP=$PATHOUT/TMP
# Limpiar posibles ejecuciones previas
#rm -r $PATHAMAS $PATHFIXCORE $PATHTREE $TMP

mkdir -p $PATHTREE
mkdir -p $PATHFIXCORE
mkdir -p $PATHCONS
mkdir -p $PATHAMAS
mkdir -p $PATHGENETREE
mkdir -p $TMP
mkdir -p $PATHALG

for file in $ALIGNMENTS/*.aln; do
# Alignment name

name=$(basename $file .aln)
logthis "Processing $name"

# Create a consensus sequence
consensus.py -n EMBOSS_001 < $file > $PATHCONS/$name.mafft.cons.fasta

# Change N with gaps in consensus
sed -i 's/n/-/g' $PATHCONS/$name.mafft.cons.fasta

# New alignment with cons in first pos
cat $PATHCONS/$name.mafft.cons.fasta $file > $PATHCONS/$name.mafft.concat.fasta

# Generate protein alignment to generate a vcf file
seqkit translate --threads 8 $PATHCONS/$name.mafft.concat.fasta  > $PATHCONS/$name.mafft.concat.aa.fasta
seqkit replace --threads 8 -p "(-)" -r 'X' -s $PATHCONS/$name.mafft.concat.aa.fasta >  $PATHCONS/$name.mafft.aaX.fasta.temp1
seqkit replace --threads 8 -p "(N)" -r '*' -s $PATHCONS/$name.mafft.aaX.fasta.temp1 >  $PATHCONS/$name.mafft.aaX.fasta
rm $PATHCONS/$name.mafft.aaX.fasta.temp1

# Generate el VCF   
snp-sites -v  $PATHCONS/$name.mafft.aaX.fasta > $PATHCONS/$name.mafft.concat.aaX.vcf

################################
#           Eval MSA           #
################################

nlines=$(grep -v "^#" $PATHCONS/$name.mafft.concat.aaX.vcf | wc -l )      
  if [[ "$nlines" -ge 3  ]]; then

    Rscript $dir/fixcore.R \
            $PATHGENETREE/$name.treefile $PATHCONS/$name.mafft.concat.aaX.vcf \
            $PATHCONS/$name.kept4.names \
            $PATHCONS/$name.removed4.names \
            $PATHALG/$name.mafft.fasta

    # Filter the alignment
    seqtk subseq $file \
        $PATHCONS/$name.kept4.names > $PATHAMAS/$name.mafft.evalmsa.fasta 

    # New alignment stats
    AMAS.py summary -i $PATHAMAS/$name.mafft.evalmsa.fasta -f fasta -d dna -s -o $PATHAMAS/$name.summary
    sed '2q;d' $PATHAMAS/$name.summary >> $PATHTREE/$PREFIX.eval4.genesummary


  else
    cp $file $PATHAMAS/$name.mafft.fasta
    AMAS.py summary -i $PATHAMAS/$name.mafft.fasta -f fasta -d dna -s -o $PATHAMAS/$name.summary
    sed '2q;d' $PATHAMAS/$name.summary >> $PATHTREE/$PREFIX.eval4.genesummary
  fi

done

# Applying trimmal
logthis "Applying trimal"
mkdir -p $PATHAMAS/trimal20missingpct
rm $PATHAMAS/*trimmal*
mv $PATHAMAS/trimal20missingpct/*fasta $PATHAMAS/

# Summary with alignments with more than 20% of missing
awk  '$6>= 20 { print $0}' $PATHTREE/$PREFIX.eval4.genesummary > $PATHAMAS/trimal20missingpct/$PREFIX.2rm.eval4.genesummary

while read -r line; do 
  gene=$(echo $line | awk '{ print $1}')
  name=$(echo $gene | sed 's/.fasta//')
  
  mv $PATHAMAS/$name* $PATHAMAS/trimal20missingpct/
  
  trimal -in $PATHAMAS/trimal20missingpct/$gene -out $PATHAMAS/$name.trimal.fasta -gt 0.8 
  sed -i 's/ .*//g' $PATHAMAS/$name.trimal.fasta
  rm $PATHAMAS/trimal20missingpct/$gene

  AMAS.py summary  -i $PATHAMAS/$name.trimal.fasta   -f fasta -d dna -s -o $PATHAMAS/$name.trimal.summary
  sed '2q;d' $PATHAMAS/$name.trimal.summary >> $PATHTREE/$PREFIX.eval4.genesummary
  
  nseqs=$(grep ">" $PATHAMAS/$name.trimal.fasta | wc -l | cut -d" " -f1)
  nlines=$(wc -l $PATHAMAS/$name.trimal.fasta | cut -d" " -f1)

  if [[ $nlines -le $nseqs ]]; then
	mv $PATHAMAS/$name.trimal.fasta $PATHAMAS/trimal20missingpct/
  fi

done <  $PATHAMAS/trimal20missingpct/$PREFIX.2rm.eval4.genesummary


awk '{ print $1 }' $PATHAMAS/trimal20missingpct/$PREFIX.2rm.eval4.genesummary > $PATHAMAS/trimal20missingpct/$PREFIX.2rm.eval4.tmp
grep -v -Ff  $PATHAMAS/trimal20missingpct/$PREFIX.2rm.eval4.tmp $PATHTREE/$PREFIX.eval4.genesummary > $PATHTREE/$PREFIX.eval4trimal.genesummary

logthis "Finished"

exit 0
