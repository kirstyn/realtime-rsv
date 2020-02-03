import csv
import sys
import argparse
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description='Trim primers.')

    parser.add_argument("--cns", action="store", type=str, dest="cns")
    parser.add_argument("--paf", action="store", type=str, dest="paf")

    parser.add_argument("--min_coverage", action="store", type=int, dest="min_coverage")

    parser.add_argument("--masked_cns", action="store", type=str, dest="masked_cns")

    return parser.parse_args()

def load_cns(cns):
    for record in SeqIO.parse(str(cns), "fasta"):
        my_record = record
    return my_record, len(my_record)

"""
to estimate coverage, for the cns- have a list [0, len(cns)] of counters that 
increments by one if there is a read mapped over that site.
Once you have that list, then iterate over the bases in the cns file 
and mask any bases as N with lower than X coverage
"""
def get_coverage(paf, cns_len):
    coverage=[]
    for i in range(cns_len):
        coverage.append(0)
    read_count = 0
    with open(paf, "r") as f:
        for l in f:
            read_count +=1
            l = l.rstrip("\n")
            tokens = l.split("\t")
            start,end = int(tokens[7]),int(tokens[8])
            map_span = range(start-1, end-1) #triple check this index thing
            for i in map_span:
                coverage[i]+=1
    print(read_count)
    return coverage

def mask_bases(coverage_list, cns_seq, min_coverage):
    masked_seq = ""
    for coverage,base in zip(coverage_list, cns_seq):
        if coverage >= min_coverage:
            masked_seq += base
        else:
            masked_seq += 'N'
    return masked_seq


if __name__ == '__main__':

    args = parse_args()

    cns, cns_len = load_cns(args.cns)

    coverage = get_coverage(args.paf,cns_len)
    masked_seq = mask_bases(coverage, cns.seq, args.min_coverage)

    with open(str(args.masked_cns), "w") as fw:
        fw.write(">{}\n{}\n".format(cns.description, masked_seq))
