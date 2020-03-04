from collections import Counter
import csv
import os
import sys
import argparse
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description='RAMPART report.')

    parser.add_argument("--annotations_path", action="store", type=str, dest="annotations_path")
    parser.add_argument("--outfile", action="store", type=str, dest="outfile")

    return parser.parse_args()


if __name__ == '__main__':

    args = parse_args()

    c_all = 0
    c_mapped = 0
    barcode_dict_counter= Counter()
    barcode_dict_mapped = Counter()

    for r,d,f in os.walk(str(args.annotations_path)):
        for fn in f:
            if fn.endswith("csv"):
                with open(r+'/'+fn,"r") as fr:
                    reader = csv.DictReader(fr)
    #                 print(reader.fieldnames)
                    for row in reader:
                        c_all +=1
                        best_reference = row["best_reference"]
                        if best_reference != "*" and best_reference != '?':
                            c_mapped +=1
                        barcode_dict_counter[row["barcode"]]+=1
                        if best_reference != "*" and best_reference != '?':
                            barcode_dict_mapped[row["barcode"]]+=1
                            
    print(f"Total reads is {c_all}\nTotal reads mapped is {c_mapped}.")
    for k in barcode_dict_counter:
        print(f"Total reads for {k} is {barcode_dict_counter[k]}\nTotal reads mapped for {k} is {barcode_dict_mapped[k]}.")
                        
                    
    with open(str(args.outfile),"w") as fw:
        fw.write("sample,total_reads,mapped_reads,unmapped_reads\n")
        fw.write(f"all,{c_all},{c_mapped},{int(c_all-c_mapped)}\n")
        for k in barcode_dict_counter:
            fw.write(f"{k},{barcode_dict_counter[k]},{barcode_dict_mapped[k]},{int(barcode_dict_counter[k]-barcode_dict_mapped[k])}\n")
            