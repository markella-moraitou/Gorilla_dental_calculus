#!/bin/bash -l

#SBATCH -A SNIC_PROJECT_ID
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -J PolyGrem_controls
#SBATCH --mail-type BEGIN
#SBATCH --mail-type END
#SBATCH --mail-type FAIL
#SBATCH --mail-user USER_EMAIL

# analysis step 1: PolyG removal -- 2nd run: for extraction, library and museum controls
# removes the polyG tails that result from the two-colour technology of the NextSeq and NovaSeq Illumina platforms 
# after concatenating different files (lane or run) for the same sample

#load modules
module load bioinfo-tools
module load fastp

echo "$SLURM_JOB_NAME"
echo $(module list)

OUTDIR=2_polyGrem

conda activate fastp #activate the environment where fastp is installed

#For DC2 library and extraction blanks (I am using only controls from the batches that contained samples from my analysis)

DATADIR=/proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_SECONDSCREEN_190219

controls=(BE108 BE208 BE109 BE209 BE110 BE210 BE111 BE211 BE112 BE212 BE213 BE114 BE214 BE115 BE215 BE116 BE216 BE117 BE217 BE118 BE218 BE119 BE219 BE319 BE121 BE221 BkL2 BL105 BL106 BL107 BL108 BL110 BL111 BL112 BS001 BS005)

for name in ${controls[@]}
do
    find $DATADIR/P1_demux_run**/**/ -name "$name*" | while read i
    do
        label=${i##/**/}
        ls -d -1 $DATADIR/P1_demux_run**/**/"$label" | xargs cat > $OUTDIR/"${label%.fastq.gz}"_concat.fastq.gz
        fastp -i $OUTDIR/"${label%.fastq.gz}"_concat.fastq.gz --trim_poly_g -A -Q -L -o $OUTDIR/"${label%.fastq.gz}"_concat_Gtrimmed.fastq.gz
        rm $OUTDIR/"${label%.fastq.gz}"_concat.fastq.gz
    done
done

#For DC2 library and extraction blanks (I am using only controls from the batches that contained samples from my analysis)

DATADIR=/proj/sllstore2017021/nobackup/ADRIAN/calculus/DC3/P1_demux_211001/

controls=(13EB2 14EB2 15EB2 17EB2 BL113 BL114 BL115 BL116 24EB1)

#forward reads
for name in ${controls[@]}
do
    find $DATADIR/** -name "$name*pair1.fastq.gz" | while read i
    do
        label=${i##/**/}
        ls -d -1 $DATADIR/**/"$label" | xargs cat > $OUTDIR/"${label%.pair1.fastq.gz}"_F_concat.fastq.gz
        fastp -i $OUTDIR/"${label%.pair1.fastq.gz}"_F_concat.fastq.gz --trim_poly_g -A -Q -L -o $OUTDIR/"${label%.pair1.fastq.gz}"_F_concat_Gtrimmed.fastq.gz
        rm $OUTDIR/"${label%.pair1.fastq.gz}"_F_concat.fastq.gz
    done
done

#reverse reads
for name in ${controls[@]}
do
    find $DATADIR/** -name "$name*pair2.fastq.gz" | while read i
    do
        label=${i##/**/}
        ls -d -1 $DATADIR/**/"$label" | xargs cat > $OUTDIR/"${label%.pair2.fastq.gz}"_R_concat.fastq.gz
        fastp -i $OUTDIR/"${label%.pair2.fastq.gz}"_R_concat.fastq.gz --trim_poly_g -A -Q -L -o $OUTDIR/"${label%.pair2.fastq.gz}"_R_concat_Gtrimmed.fastq.gz
        rm $OUTDIR/"${label%.pair2.fastq.gz}"_R_concat.fastq.gz
    done
done

#For Jena library and extraction blanks

DATADIR=/proj/sllstore2017021/nobackup/GORILLA_METAGENOMES/blanks

#forward reads
find $DATADIR -name "*_S0_L000_R1_000.trimmed.fastq.gz" | while read i
do
    label=${i##/**/}
    fastp -i "$i" --trim_poly_g -A -Q -L -o $OUTDIR/"${label%_S0_L000_R1_000.trimmed.fastq.gz}"_F_concat_Gtrimmed.fastq.gz
done

#reverse reads
find $DATADIR -name "*_S0_L000_R2_000.trimmed.fastq.gz" | while read i
do
    label=${i##/**/}
    fastp -i "$i" --trim_poly_g -A -Q -L -o $OUTDIR/"${label%_S0_L000_R2_000.trimmed.fastq.gz}"_R_concat_Gtrimmed.fastq.gz
done


#Museum controls and DC1 blanks were sequenced in HiSeq so they will be skipped in this step

#Rename 24EB1 to make it match DC2 blanks
ls 24EB1* | while read i
do
    mv "$i" BE124"${i#24EB1}"
done

conda deactivate

