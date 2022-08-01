#Adapted from SNIC_PROJECT/SAM/Bear_Lineage_Project/scripts/DistanceToReferences.py

#Import the modules to be used
import os
import pysam #this is a module for parsing alignment files, will need to be installed
import argparse

############ Defining the functions we will use

#This is the function that open and transforms the fasta file to a dictionary in its simplest use (it also subsets by position and sample if wanted, that's why its long and messy.).
def FastaFetch(fastafile, sequences=None, start=None, stop=None):
    f = pysam.FastaFile(fastafile)
    seqDict = []
    if sequences is not None:
        if os.path.isfile(sequences):
            seqList = []
            with open(sequences) as s:
                for line in s.readlines():
                    if not line == "":
                        seqList.append(line.rstrip().lstrip(">")) 
        else:
            seqList = []
            for s in sequences.split(","):
                if not s == "":
                    seqList.append(s.rstrip().lstrup(">"))
        if start and stop:
            for seq in seqList:
                seqDict.append({'header':seq, 'sequence':f.fetch(seq, start, stop)})
        elif start:
            for seq in seqList:
                seqDict.append({'header':seq, 'sequence':f.fetch(seq, start)})
        else:
            for seq in seqList:
                seqDict.append({'header':seq, 'sequence':f.fetch(seq)})
    else:
        if start and stop:
            seqDict = []
            for seq in f.references:
                seqDict.append({'header':seq, 'sequence':f.fetch(seq, start, stop)})
        elif start:
            for seq in f.references:
                seqDict.append({'header':seq, 'sequence':f.fetch(seq, start)})
        else:
            for seq in f.references:
                seqDict.append({'header':seq, 'sequence':f.fetch(seq)})
    return seqDict

#This function checks pop1's genotypes and when that's fixed, compares it to pop2's genotype and outputs the position and bases if different
def CountDiagnosticPositions(seqDict, reference, outfile):
    #First line should hold the popnames
    with open(reference) as f:
        header = f.readline().split("\t")
        pop1 = header[1]
        pop2 = header[2].rstrip()
        #and line 1- should have the sites
        f.seek(0)
        sites = f.read().split("\n")[1:]
        distanceDict = []
        for seq in seqDict:
            #count number of similar sites per sample, start at 0 for each sample
            pop1_similarity = 0
            pop2_similarity = 0
            other = 0
            missing_data = 0
            for site in sites:
                if site.split("\t")[0] == "": #if position is a blank line, break and start with next sample
                    break
                pos = int(site.split("\t")[0]) - 1 #zero based position in python, so if input is not one-based this needs changing (it is one-based when using FindDiagnosticSites.py output)
                pop1_allele = site.split("\t")[1]
                pop2_allele = site.split("\t")[2]
                #count the similarities: for now it's still case sensitive! 
                if seq['sequence'][pos] == pop1_allele:
                    pop1_similarity += 1
                elif seq['sequence'][pos] == pop2_allele:
                    pop2_similarity += 1
                elif seq['sequence'][pos] in ['N', 'n', '-', '?']:
                    missing_data += 1
                else:
                    other += 1
            #add the samples stats to dictionary
            distanceDict.append({'sample':seq['header'], pop1:pop1_similarity,
            pop2:pop2_similarity, 'other':other, 'missing_data':missing_data})
    #then we write those to a tsv outfile:
    with open(outfile, "w") as of:
        of.write("Sample\t" + pop1 + "\t" + pop2 + "\t" + "other" + "\t" + "missing" + "\n")
        for sample in distanceDict:
            of.write(sample['sample'] + "\t" + str(sample[pop1]) + "\t" + str(sample[pop2]) + "\t" + str(sample['other']) + "\t" + str(sample['missing_data']) + "\n")
        of.close()

###Input arguments
parser = argparse.ArgumentParser(description="Count diagnostic sites per sample. Diagnostic sites reference file given as input.")
parser.add_argument('-i', '--input', type=str, help='Input fasta alignment.', required=True)
parser.add_argument('--reference', type=str, help="File with reference's diagnostic sites, as output by the script FindDiagnosticSites.py.", required=True)
parser.add_argument('-o', '--output', type=str, help="Path to output file containing the diagnostic counts per sample", required=True)

args = parser.parse_args()

## First handle the input fasta to a dictionary
seqDict = FastaFetch(args.input)

#then count all diagnostic sites per sample and output to tsv
CountDiagnosticPositions(seqDict, args.reference, args.output)


