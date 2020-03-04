import argparse
from Bio import SeqIO
import collections
import sys

def parse_args():
    parser = argparse.ArgumentParser(description='Parse mappings, add to headings and create report.')

    parser.add_argument("--csv", action="store", type=str, dest="in_csv")
    parser.add_argument("--amplicons", action="store", type=str, dest="amplicons")
    parser.add_argument("--out_csv", action="store", type=str, dest="out_csv")

    return parser.parse_args()

def make_amp_dict(amplicons):
    amp_dict = collections.defaultdict(list)
    with open(str(amplicons), "r") as f:
        for l in f:
            l=l.rstrip("\n")
            tokens= l.split(",")
            accession,reference,amplicon= tokens[0].split('|')
            coords=tokens[1].split(':')
            if amplicon in ["Amp1","Amp2","Amp3"]:
                amp_dict["Amp123"].append(int(coords[0]))
                amp_dict["Amp123"].append(int(coords[1]))
            if amplicon in ["Amp4"]:
                amp_dict["Amp4"].append(int(coords[0]))
                amp_dict["Amp4"].append(int(coords[1]))
            if amplicon in ["Amp5"]:
                amp_dict["Amp5"].append(int(coords[0]))
                amp_dict["Amp5"].append(int(coords[1]))
    amp_coords = {}
    for amp in amp_dict:
        s= sorted(amp_dict[amp])
        start,end = s[0],s[-1]
        amp_coords[amp] = (start,end)
    return amp_coords

def check_overlap(coords1,coords2):
    list1 = list(range(coords1[0],coords1[1]))
    list2 = list(range(coords2[0],coords2[1]))
    overlap = set(list1).intersection(list2)
    if overlap:
        return True, len(overlap)
    else:
        return False, 0 

def parse_line(line):

    values = {}

    tokens = line.rstrip('\n').split(',')
    values["read_name"], values["read_len"] = tokens[:2]
    values["ref_hit"], values["ref_len"], values["coord_start"], values["coord_end"], values["matches"], values["aln_block_len"] = tokens[4:10]
    return values

def parse_csv(in_csv, amp_dict, out_csv):
    counts = {
        "unmapped": 0,
        "total": 0
    }
    
    with open(str(in_csv),"r") as fr:
        for line in fr:
            line = line.rstrip('\n')
            if line.startswith("read_name"):
                csv_report.write(f"{line},amplicon\n")
            else:
                counts["total"]+=1
                mapping = parse_line(line)
                if mapping["ref_hit"] in ['*','?','none']:
                    counts["unmapped"]+=1
                    out_csv.write(line+",nan\n")
                else:
                    overlap_list = []
                    for i in amp_dict:
                        overlap, length = check_overlap((int(mapping["coord_start"]),int(mapping["coord_end"])),amp_dict[i])
                        if overlap:
                            overlap_list.append((i, length))
                    amplicon = sorted(overlap_list, key = lambda x : x[1], reverse=True)[0][0]
                    out_csv.write(line+f",{amplicon}\n")

if __name__ == '__main__':

    args = parse_args()

    with open(str(args.out_csv), "w") as csv_report:
        amp_dict = make_amp_dict(args.amplicons)
        parse_csv(args.in_csv, amp_dict, csv_report)