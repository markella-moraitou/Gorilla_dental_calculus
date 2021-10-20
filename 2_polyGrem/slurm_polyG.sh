#!/bin/bash -l

#SBATCH -A snic2020-5-528
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 6:00:00
#SBATCH -J PolyGrem
#SBATCH --mail-type BEGIN
#SBATCH --mail-type END
#SBATCH --mail-type FAIL
#SBATCH --mail-user Markella.Moraitou.0437@student.uu.se

# analysis step 1: PolyG removal
# removes the polyG tails that result from the two-colour technology of the NextSeq and NovaSeq Illumina platforms 
# after concatenating different files (lane or run) for the same sample

#load modules
module load bioinfo-tools
module load fastp

OUTDIR=/proj/sllstore2017021/nobackup/MARKELLA/2_polyGrem

conda activate fastp #activate the environment where fastp is installed

#For DC2 samples

DATADIR=/proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_SECONDSCREEN_190219/

find $DATADIR/P1_demux_run**/**/ -name "G**.fastq.gz" | while read i
do
    name=${i##/**/} #remove path
    echo $name
done | sort | uniq | while read n #iterate through list of unique names
do #concatenated them when they have the same name, pipe through fastp to trim polyG and save in OUTDIR
    ls -d -1 $DATADIR/P1_demux_run**/**/$n | xargs cat > $OUTDIR/${n%.fastq.gz}_concat.fastq.gz
    fastp -i $OUTDIR/${n%.fastq.gz}_concat.fastq.gz --trim_poly_g -A -Q -L -o $OUTDIR/${n%.fastq.gz}_concat_Gtrimmed.fastq.gz
    rm $OUTDIR/${n%.fastq.gz}_concat.fastq.gz
done

#For DC3 samples

DATADIR=/proj/sllstore2017021/nobackup/ADRIAN/calculus/DC3/P1_demux_211001/

#forward reads
find $DATADIR/**/ -name "G**pair1**.fastq.gz" | while read i
do
    oldname=${i##/**/} #remove path
    echo $oldname
done | sort | uniq | while read n #iterate through list of unique names
do #concatenated them when they have the same name, pipe through fastp to trim polyG and save in OUTDIR
    #Fix name to match DC2 names
    number=${n#G}
    number=${number%%.*}
    #Change the name to match DC2
    if [[ ${number%_*} -lt 10 ]]
    then
        name=G000${number}_F.fastq.gz
    else
        name=G00${number}_F.fastq.gz
    fi
    ls -d -1 $DATADIR/**/$n | xargs cat > $OUTDIR/${name%.fastq.gz}_concat.fastq.gz
    fastp -i $OUTDIR/${name%.fastq.gz}_concat.fastq.gz --trim_poly_g -A -Q -L -o $OUTDIR/${name%.fastq.gz}_concat_Gtrimmed.fastq.gz
    rm $OUTDIR/${name%.fastq.gz}_concat.fastq.gz
done

#reverse reads
find $DATADIR/**/ -name "G**pair2**.fastq.gz" | while read i
do
    oldname=${i##/**/} #remove path
    echo $oldname
done | sort | uniq | while read n #iterate through list of unique names
do #concatenated them when they have the same name, pipe through fastp to trim polyG and save in OUTDIR
    #Fix name to match DC2 names
    number=${n#G}
    number=${number%%.*}
    #Change the name to match DC2
    if [[ ${number%_*} -lt 10 ]]
    then
        name=G000${number}_R.fastq.gz
    else
        name=G00${number}_R.fastq.gz
    fi
    ls -d -1 $DATADIR/**/$n | xargs cat > $OUTDIR/${name%.fastq.gz}_concat.fastq.gz
    fastp -i $OUTDIR/${name%.fastq.gz}_concat.fastq.gz --trim_poly_g -A -Q -L -o $OUTDIR/${name%.fastq.gz}_concat_Gtrimmed.fastq.gz
    rm $OUTDIR/${name%.fastq.gz}_concat.fastq.gz
done

#For Jena samples

DATADIR=/proj/sllstore2017021/nobackup/GORILLA_METAGENOMES/pairedend/paired_fastq #incomplete path: there are directories paired_fastq2 and paired_fastq3

#forward reads
find ${DATADIR}[23]/* -name "*R1_000.trimmed.fastq.gz" | while read i
do
    name=${i##/**/} #remove path
    name=${name:0:6} #keep only extraction ID
    echo $name
done | sort | uniq | while read n #iterate through list of unique names
do #concatenated them when they have the same name, pipe through fastp to trim polyG and save in OUTDIR (excluding UDG-treated samples)
    ls -d -1 ${DATADIR}[23]/*/${n}.A0101**R1_000.trimmed.fastq.gz | xargs cat > $OUTDIR/${n}_F_concat.fastq.gz
    fastp -i $OUTDIR/${n}_F_concat.fastq.gz --trim_poly_g -A -Q -L -o $OUTDIR/${n}_F_concat_Gtrimmed.fastq.gz
    rm $OUTDIR/${n}_F_concat.fastq.gz
done

#reverse reads
find ${DATADIR}[23]/* -name "*R2_000.trimmed.fastq.gz" | while read i
do
    name=${i##/**/} #remove path
    name=${name:0:6} #keep only extraction ID
    echo $name
done | sort | uniq | while read n #iterate through list of unique names
do #concatenated them when they have the same name, pipe through fastp to trim polyG and save in OUTDIR (excluding UDG-treated samples)
    ls -d -1 ${DATADIR}[23]/*/${n}.A0101**R2_000.trimmed.fastq.gz | xargs cat > $OUTDIR/${n}_R_concat.fastq.gz
    fastp -i $OUTDIR/${n}_R_concat.fastq.gz --trim_poly_g -A -Q -L -o $OUTDIR/${n}_R_concat_Gtrimmed.fastq.gz
    rm $OUTDIR/${n}_R_concat.fastq.gz
done
