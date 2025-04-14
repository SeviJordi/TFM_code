#!/bin/bash
#SBATCH --job=snps
#SBATCH --output=snps_%j.out
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=50
#SBATCH --time=05-00:00:00
#SBATCH --mem=5G
#SBATCH --qos=medium
#SBATCH  -N 3

# Jordi Sevilla Fortuny
# Script to calculate SNPs from between homologous sequences of different approaches
# Usage: JSF_TFM_06_snps.sh [gene names] [subset path] [ST] [output dir]
# Input: gene names file, subset path, ST, output directory
# Output: VCF files with SNPs for each gene family

# Check arguments
if [ "$#" -ne 4 ]; then
  echo "Incorrect number of arguments"
  echo "Usage: $0 [gene names] [subset path] [ST] [output dir]"
  exit 1
fi

gene_names=$1 # Input file with correspondance between gene families of different approaches
subset_path=$2 # Path to were the results of subset pangenomes are
ST=$3
outdir=$4
tmpdir="tmp/$ST"

mkdir -p "$tmpdir"

# Create output
mkdir -p $outdir

# Analyze
mkdir -p $outdir/$ST

calculate_div(){  # line #subset_path #ST #outdir $tmpdir
	ST=$3
	subset_path=$2
	line=$1
	outdir=$4
	
	while read name panacota panaroo panaroodefault roarydefault; do
		
		new_panacota=$(echo $panacota | sed 's/'$ST'/'$ST'-/g') # Fix panacota names
	        panacota=$new_panacota

        # Select the snippy alignment depending on the avalible pangenomes for the gene family
		 if [[ $panacota != "NA" ]]; then
                        snippy_al=$subset_path/${ST}_panacota/snippy_genes/$panacota.aln
                elif [[ $panaroo != "NA" ]]; then
                        snippy_al=$subset_path/${ST}_panaroo/snippy_genes/$panaroo.aln
                elif [[ $panaroodefault != "NA" ]]; then
                        snippy_al=$subset_path/${ST}_panaroodefault/snippy_genes/$panaroodefault.aln
                elif [[ $roarydefault != "NA" ]]; then
                        snippy_al=$subset_path/${ST}_roarydefault/snippy_genes/$roarydefault.aln
                fi

			for seq in $(grep ">" $snippy_al | sed 's/>//g'); do
				
				mafft --thread 2 --adjustdirection <(
                                                cat \
							                            <(echo $seq | seqtk subseq $snippy_al - | sed 's/>.*/>'$seq'_fakeref/g') \
		                                                <(echo $seq | seqtk subseq $snippy_al - | sed 's/>.*/>'$seq'_snippy/g') \
                                                        <(echo $seq | seqtk subseq $subset_path/${ST}_panacota/raw_genes/$panacota.aln - | sed 's/>.*/>'$seq'_pc/g') \
                                                        <(echo $seq | seqtk subseq $subset_path/${ST}_panacota/clipkit_genes/$panacota.aln - | sed 's/>.*/>'$seq'_pc_c/g') \
                                                        <(echo $seq | seqtk subseq $subset_path/${ST}_panacota/trimal_genes/$panacota.aln - | sed 's/>.*/>'$seq'_pc_t/g') \
                                                        <(echo $seq | seqtk subseq $subset_path/${ST}_panacota/treeshrink_genes/$panacota.aln - | sed 's/>.*/>'$seq'_pc_ts/g') \
                                                        <(echo $seq | seqtk subseq $subset_path/${ST}_panacota/hmmcleaner_genes/$panacota.aln - | sed 's/>.*/>'$seq'_pc_h/g') \
                                                        <(echo $seq | seqtk subseq $subset_path/${ST}_panacota/fixcore_genes/$panacota.aln - | sed 's/>.*/>'$seq'_pc_f/g') \
                                                        <(echo $seq | seqtk subseq $subset_path/${ST}_panaroo/raw_genes/$panaroo.aln - | sed 's/>.*/>'$seq'_pr/g') \
                                                        <(echo $seq | seqtk subseq $subset_path/${ST}_panaroo/clipkit_genes/$panaroo.aln - | sed 's/>.*/>'$seq'_pr_c/g') \
                                                        <(echo $seq | seqtk subseq $subset_path/${ST}_panaroo/trimal_genes/$panaroo.aln - | sed 's/>.*/>'$seq'_pr_t/g') \
                                                        <(echo $seq | seqtk subseq $subset_path/${ST}_panaroo/treeshrink_genes/$panaroo.aln - | sed 's/>.*/>'$seq'_pr_ts/g') \
                                                        <(echo $seq | seqtk subseq $subset_path/${ST}_panaroo/hmmcleaner_genes/$panaroo.aln - | sed 's/>.*/>'$seq'_pr_h/g') \
                                                        <(echo $seq | seqtk subseq $subset_path/${ST}_panaroo/fixcore_genes/$panaroo.aln - | sed 's/>.*/>'$seq'_pr_f/g') \
                                                        <(echo $seq | seqtk subseq $subset_path/${ST}_panaroodefault/raw_genes/$panaroodefault.aln - | sed 's/>.*/>'$seq'_prf/g') \
                                                        <(echo $seq | seqtk subseq $subset_path/${ST}_panaroodefault/clipkit_genes/$panaroodefault.aln - | sed 's/>.*/>'$seq'_prf_c/g') \
                                                        <(echo $seq | seqtk subseq $subset_path/${ST}_panaroodefault/trimal_genes/$panaroodefault.aln - | sed 's/>.*/>'$seq'_prf_t/g') \
                                                        <(echo $seq | seqtk subseq $subset_path/${ST}_panaroodefault/treeshrink_genes/$panaroodefault.aln - | sed 's/>.*/>'$seq'_prf_ts/g') \
                                                        <(echo $seq | seqtk subseq $subset_path/${ST}_panaroodefault/hmmcleaner_genes/$panaroodefault.aln - | sed 's/>.*/>'$seq'_prf_h/g') \
                                                        <(echo $seq | seqtk subseq $subset_path/${ST}_panaroodefault/fixcore_genes/$panaroodefault.aln - | sed 's/>.*/>'$seq'_prf_f/g') \
                                                        <(echo $seq | seqtk subseq $subset_path/${ST}_roarydefault/raw_genes/$roarydefault.aln - | sed 's/>.*/>'$seq'_rf/g') \
                                                        <(echo $seq | seqtk subseq $subset_path/${ST}_roarydefault/clipkit_genes/$roarydefault.aln - | sed 's/>.*/>'$seq'_rf_c/g') \
                                                        <(echo $seq | seqtk subseq $subset_path/${ST}_roarydefault/trimal_genes/$roarydefault.aln - | sed 's/>.*/>'$seq'_rf_t/g') \
                                                        <(echo $seq | seqtk subseq $subset_path/${ST}_roarydefault/treeshrink_genes/$roarydefault.aln - | sed 's/>.*/>'$seq'_rf_ts/g') \
                                                        <(echo $seq | seqtk subseq $subset_path/${ST}_roarydefault/hmmcleaner_genes/$roarydefault.aln - | sed 's/>.*/>'$seq'_rf_h/g') \
                                                        <(echo $seq | seqtk subseq $subset_path/${ST}_roarydefault/fixcore_genes/$roarydefault.aln - | sed 's/>.*/>'$seq'_rf_f/g') \
                                        ) > $5/$name@$seq.aln

				snp-sites -v $5/$name@$seq.aln  > $4/$name@$seq.vcf
			if [[ ! -s $4/$name@$seq.vcf ]]; then
				rm $4/$name@$seq.vcf
			fi
 			done	
		
	done < <(echo $line)

	
	    	

}

export -f calculate_div

parallel -j 150 -a $gene_names calculate_div {}  $subset_path $ST $outdir/$ST $tmpdir

rm -rf $tmpdir

exit 0