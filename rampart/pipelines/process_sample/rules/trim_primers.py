import csv
import sys
import argparse
from Bio.SeqIO.QualityIO import FastqGeneralIterator
# at the moment just trims off 28 bases from the start and end of a read. 
# Need to adapt to trimming more specifically


def parse_args():
    parser = argparse.ArgumentParser(description='Trim primers.')
    parser.add_argument("--reads", action="store", type=str, dest="reads")

    parser.add_argument("--output_reads", action="store", type=str, dest="output_reads")

    return parser.parse_args()

def file_writer(reads,output_reads,trim):
    trim = int(trim)
    with open(output_reads, "w") as handle:
        for title, seq, qual in FastqGeneralIterator(open(reads)):
            handle.write(f"@{title}\n{seq[trim:-trim]}\n+\n{qual[trim:-trim]}\n")


if __name__ == '__main__':

    args = parse_args()

    file_writer(args.reads,args.output_reads,28)
