#Adapted from /proj/sllstore2017021/nobackup/SAM/Bear_Lineage_Project/scripts/FindDiagnosticSites_edited.py

###This script will identify positions in an alignment that are fixed for differential alleles between two populations
###Takes a fasta alignment as input, and a file assigning sequences to one of two populations
###The popfile should have the sample/header name in the first column and the population assignment in the second.
###Samples that doesn't have a population assignment will be ignored in this analysis

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
                    seqList.append(s.rstrip().lstrip(">"))
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
def FindDiagnosticPositions(seqDict, popfile, outfile):
    with open(popfile) as pf:
        popmap = pf.read().split("\n")
        popdict = []
        for sample in popmap:
            if len(sample.split("\t")) == 2:
                popdict.append({'header':sample.split("\t")[0], 'pop':sample.split("\t")[1]})
    
    #a dictionary with pop as key and samples assigned as values
    #first get the names of the populations
    pop1 = popdict[0]['pop']
    for pop in popdict:
        if not pop['pop'] == pop1:
            pop2 = pop['pop']
            break
    #then put in the samples:
    popmap = {pop1: [], pop2: []}
    for i in popdict:
        if i['pop'] == pop1:
            popmap[pop1].append(i['header'].lstrip(">"))
        elif i['pop'] == pop2:
            popmap[pop2].append(i['header'].lstrip(">"))
    print("assigning the following samples as references: ")
    print(popmap)
    #go through all the bases in the alignment and check them
    diagnosticSites = []
    for base in range(0, len(seqDict[0]['sequence'])):
        pop1_alleles = []
        pop2_alleles = []
        for seq in seqDict:
            if seq['header'] in popmap[pop1]:
                if not seq['sequence'][base] in ["N","n","-","?"]:
                    pop1_alleles.append(seq['sequence'][base])
            elif seq['header'] in popmap[pop2]:
                if not seq['sequence'][base] in ["N","n","-","?"]:
                    pop2_alleles.append(seq['sequence'][base])
        if len(set(pop1_alleles)) == 1:
            if len(set(pop2_alleles)) == 1:
                if not pop1_alleles[0] == pop2_alleles[0]:
                    diagnosticSites.append((base, {pop1:pop1_alleles[0], pop2:pop2_alleles[0]}))
    if len(diagnosticSites) == 0:
        print("No fixed differences were found")
    #then we write those to a tsv outfile:
    with open(outfile, "w") as of:
        of.write("Position\t" + pop1 + "\t" + pop2 + "\n")
        for site in diagnosticSites:
            of.write(str(site[0] + 1) + "\t" + site[1][pop1] + "\t" + site[1][pop2] + "\n")
        of.close()


#Input arguments
parser = argparse.ArgumentParser(description="Identify sites that are fixed for differential alleles between two predefined populations")
parser.add_argument('-i', '--input', type=str, help='Input fasta alignment.', required=True)
parser.add_argument('--popfile', type=str, help="File with two columns separated by a tab. First column should have sample name and second the population it belongs to. Second column can be empty and then sample will be ignored.")
parser.add_argument('-o', '--output', type=str, help="Path to output file containing the fixed differences.")

args = parser.parse_args()


### Now we run the functions we've defined above

## First handle the input fasta to a dictionary
seqDict = FastaFetch(args.input)
print("Found " + str(len(seqDict)) + " sequences in the alignmant.")
#then run this seqdict through the next function along with the popfile and output from input arguments
FindDiagnosticPositions(seqDict, args.popfile, args.output)

print("First sequence length: " + str(len(seqDict[0]['sequence'])))
