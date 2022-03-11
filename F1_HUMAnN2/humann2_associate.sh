#!/bin/bash
# post humann analyses

module load bioinfo-tools HUMAnN/3.0 R_packages

DATADIR=F1_HUMAnN2
OUTDIR=F1_HUMAnN2/association_test
PLOTDIR=F1_HUMAnN2/plots
threads=$SLURM_NTASKS

# for association test and plotting, humann3 has better functions
# so purge modules and load again

# convert tab seperated to biom format
biom convert -i $DATADIR/all_GOterms_cpm_renamed.tsv -o $OUTDIR/all_GOterms_cpm_renamed.biom --table-type "Gene table" --to-hdf5

# strip header
head -n 1 $DATADIR/all_GOterms_cpm_renamed.tsv | tr "\t" "\n" > $OUTDIR/all_GOterms_header.tsv

cd $OUTDIR
# process header and add metadata
Rscript --vanilla addMetadata_adrian.R all_GOterms_header.tsv all_GOterms_header_wfeature.tsv

# insert on second line, with appropriate first column
sed "s/UNMAPPED/ $(<$OUTDIR/all_GOterms_header_wfeature.tsv)\nUNMAPPED/" $DATADIR/all_GOterms_cpm_renamed.tsv > all_GOterms_cpm_renamed_metadata.tsv

# run association test
# output format:
# feature, level means, p-value, q-value
humann_associate --input all_GOterms_cpm_renamed_metadata.tsv --last-metadatum " plot.label" --focal-metadatum " plot.label" --focal-type categorical --output GO_terms_stats.txt

# take significant results from association test and make plots:
Rscript --vanilla sortSigOutput_adrian.R GO_terms_stats.txt sig_GO_terms.txt

grep -v "^feature" sig_GO_terms.txt | awk '{print $2}' | head -10 | while read line
do
    # make barplot
    humann_barplot --input all_GOterms_cpm_renamed_metadata.tsv --focal-feature $line --focal-metadatum " Spec.subspecies" -s sum metadata --output $PLOTDIR/GO_terms_${line}_bar.png
done
