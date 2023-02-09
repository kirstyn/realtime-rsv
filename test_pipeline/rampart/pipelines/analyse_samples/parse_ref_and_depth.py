import argparse
from Bio import SeqIO
from collections import OrderedDict
from collections import Counter
from collections import defaultdict
import pandas as pd

import sys
import csv


def parse_args():
    parser = argparse.ArgumentParser(description='Parse mappings, add to headings and create report.')

    parser.add_argument("--csv", action="store", type=str, dest="csv")
    parser.add_argument("--reads", action="store", type=str, dest="reads")
    parser.add_argument("--references", action="store", type=str, dest="references")
    parser.add_argument("--out_counts", action="store", type=str, dest="out_counts")
    parser.add_argument("--sample", action="store", type=str, dest="sample")
    parser.add_argument("--min_reads", action="store", type=int, dest="min_reads")
    parser.add_argument("--min_pcent", action="store", type=float, dest="min_pcent")

    parser.add_argument("--output_path", action="store", type=str, dest="output_path")

    return parser.parse_args()

def make_ref_dict(references):
    refs = {}
    for record in SeqIO.parse(references,"fasta"):
        refs[record.id]=record.seq
    return refs

def count_and_return_analysis_dict(report,csv_out,sample):
    counts = OrderedDict()
    counts["RABV_PER"]=0
    counts["RABV_BRA"]=0
    counts["RABV_PHL"]=0
    counts["RABV_TZA"]=0
    counts["*"]=0
    counts["?"]=0

    detail_dict= {
        "RABV_PER": Counter(),
        "RABV_BRA": Counter(),
        "RABV_PHL": Counter(),
        "RABV_TZA": Counter(),
        "*": Counter(),
        "?": Counter()
    }

    total = 0
    with open(str(report),"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            total +=1
            detail_dict[row["display_name"]][row["best_reference"]]+=1
            counts[row["display_name"]]+=1

    count_str = f"{sample},"
    for i in counts:
        count_str += f"{counts[i]},"
    count_str = count_str.rstrip(',') + '\n'
    csv_out.write(count_str)

    analysis_dict = {}
    for key in counts:
        if key not in ['*',"?"]:
            if counts[key] > args.min_reads and 100*(counts[key]/total)> args.min_pcent:
            
                top = detail_dict[key].most_common()[0]
                analysis_dict[key] = (top[0],f">{key} best_reference={top[0]} num_reads={counts[key]}")

    read_dict = defaultdict(list)

    with open(str(report),"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["display_name"] in analysis_dict:
                read_dict[row["display_name"]].append(row["read_name"])
    print(','.join(analysis_dict.keys()))

    return analysis_dict, read_dict
    

if __name__ == '__main__':

    args = parse_args()

    ref_dict = make_ref_dict(str(args.references))

    csv_report = open(str(args.out_counts), "w")
    analysis_dict,read_dict = count_and_return_analysis_dict(args.csv, csv_report, args.sample)

    for ref in analysis_dict:
        ref_file = args.output_path + "/" + ref + ".fasta"
        best_ref = analysis_dict[ref][0]
        header = analysis_dict[ref][1]
        with open(ref_file,"w") as fw:
            fw.write(f"{header}\n{ref_dict[best_ref]}\n")
        
        read_file = args.output_path + "/" + ref + ".fastq"
        with open(read_file,"w") as fw:
            records = []
            for record in SeqIO.parse(args.reads,"fastq"):
                if record.id in read_dict[ref]:
                    records.append(record)

            SeqIO.write(records, fw, "fastq")

    csv_report.close()

