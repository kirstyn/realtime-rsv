import csv
import sys
import argparse
from Bio import AlignIO

def parse_args():
    parser = argparse.ArgumentParser(description='Trim primers.')

    parser.add_argument("-o", action="store", type=str, dest="o")
    parser.add_argument("-i", action="store", type=str, dest="i")

    return parser.parse_args()


def merge_seqs(aln):
    alignment = AlignIO.read(aln, "fasta")
    cns_string = ""
    for i in range(len(alignment[0])):
        col = alignment[:,i]
        if '-' in col:
            base = [i for i in col if i!= '-'][0]
        elif len(set(col)) == 1:
            base = col[0]
        elif len(set(col)) > 1:
            if 'n' in col.lower():
                base = [i for i in col if i!= 'n'][0]
            else:
                base = 'n'
        cns_string += base
    header =''
    accession = ''
    barcode = ''
    for i in alignment:
        info = i.description.split(" ")
        name = info[0].split("_")
        barcode = f"{name[0]} "
        genotype = f"{name[1]}={name[2]} "
        header += genotype
        accession += info[1].split("=")[1] + ':'
    accession = accession.rstrip(":")
    header= f"{barcode}{header}accession={accession} length={ len(cns_string)}"
    return header, cns_string

if __name__ == '__main__':

    args = parse_args()

    header, seq = merge_seqs(str(args.i))

    with open(str(args.o), "w") as fw:
        fw.write(">{}\n{}\n".format(header, seq))
