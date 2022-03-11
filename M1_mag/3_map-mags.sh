#!/bin/bash -l
#SBATCH -A SNIC_PROJECT_ID
#SBATCH -p core
#SBATCH --ntasks 8 --mem-per-cpu=2G
#SBATCH --array=1-46%20
#SBATCH -t 00:15:00
#SBATCH -J magmap
#SBATCH --mail-type=FAIL
#SBATCH -o slurm-%A_%a.out

# align MAGs
module load bioinfo-tools samtools bbmap vcftools bcftools

threads=$SLURM_NTASKS

# INPUTS

BASEDIR=/home/adrianf/project-folder/nobackup/ADRIAN/calculus/Gorilla/mag/metaWRAP/full-decontam-set
echo "$BASEDIR"

BIN_DIR=${BASEDIR}/MAG_binsrefined_metaWRAP/metawrap_50_10_bins
echo "$BIN_DIR"

REF=$BIN_DIR/bin.18.fa
echo "$REF"

INTERLEAVED=/home/adrianf/project-folder/nobackup/MARKELLA/RD2_mapping/decontam_files
echo "$INTERLEAVED"

OUTDIR=$BASEDIR/phylophlan/Lactobacillus
echo "$OUTDIR"
mkdir -p "$OUTDIR"

i=$(find "$INTERLEAVED"/*_m_decontam.fastq.gz | grep -v -E "BL|BE|BS|ERR" | awk -F "/" '{print $NF}' | sed -n "$SLURM_ARRAY_TASK_ID"p)

#view all mapped reads | sort ready for binning
# convert, sort, and turn into a fasta
bbmap.sh ref="$REF" in="$INTERLEAVED"/"$i" out="$OUTDIR"/"${i%_m_decontam.fastq.gz}"_mappedcontigs.sam maxindel=80 pigz=t unpigz=t

samtools view -Sb -F 4 -@ "$threads" "$OUTDIR"/"${i%_m_decontam.fastq.gz}"_mappedcontigs.sam | samtools sort -@ "$threads" - > "$OUTDIR"/"${i%_m_decontam.fastq.gz}"_mappedcontigs.sorted.bam

samtools mpileup -uf "$REF" "$OUTDIR"/"${i%_m_decontam.fastq.gz}"_mappedcontigs.sorted.bam | bcftools call -mv -Oz -o "$OUTDIR"/"${i%_m_decontam.fastq.gz}"_mappedcontigs.sorted.vcf.gz

tabix "$OUTDIR"/"${i%_m_decontam.fastq.gz}"_mappedcontigs.sorted.vcf.gz 

cat "$REF" | bcftools consensus "$OUTDIR"/"${i%_m_decontam.fastq.gz}"_mappedcontigs.sorted.vcf.gz  > "$OUTDIR"/"${i%_m_decontam.fastq.gz}"_mappedcontigs.sorted.fa

cat "$REF" | vcf-consensus "$OUTDIR"/"${i%_m_decontam.fastq.gz}"_mappedcontigs.sorted.vcf.gz  > "$OUTDIR"/"${i%_m_decontam.fastq.gz}"_mappedcontigs.sorted.fa

if [[ -f "$OUTDIR/${i%_m_decontam.fastq.gz}_mappedcontigs.sorted.fa" && -s "$OUTDIR/${i%_m_decontam.fastq.gz}_mappedcontigs.sorted.fa" ]]; then
      rm -f "$OUTDIR"/"${i%_m_decontam.fastq.gz}"_mappedcontigs.sam "$OUTDIR"/"${i%_m_decontam.fastq.gz}"_mappedcontigs.sorted.bam "$OUTDIR"/"${i%_m_decontam.fastq.gz}"_mappedcontigs.sorted.vcf*
    else
      echo "consensus sequence generated"
    fi
