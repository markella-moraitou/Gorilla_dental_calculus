#!/usr/bin/python
from sys import argv
from itertools import izip, islice
import gzip

def trim_fwd_barcodes(FASTQ_FILE,basename):
    outputfile = open(str(basename)+"_bctrim.fastq", "w")
    fastq_dict = {}

    with gzip.open(FASTQ_FILE, 'r') as readfile:
        while True:
            fastq_read_F = list(islice(readfile, 4))
            if not fastq_read_F:
                break
            if fastq_read_F[1].strip() in fastq_dict:
                continue
            else:
                fastq_dict[fastq_read_F[1].strip()] = ""
            outputfile.write(fastq_read_F[0] + fastq_read_F[1][7:] + fastq_read_F[2] + fastq_read_F[3][7:])

        outputfile.close()

if __name__ == "__main__":
    FASTQ_FILE = argv[1]
    basename = argv[2]
    trim_fwd_barcodes(FASTQ_FILE,basename)
