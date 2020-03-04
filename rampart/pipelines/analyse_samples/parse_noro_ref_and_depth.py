import argparse
from Bio import SeqIO
import collections
from collections import defaultdict
import csv
import sys

def parse_args():
    parser = argparse.ArgumentParser(description='Assess reference depth per amplicon/ orf and bin appropriately.')

    parser.add_argument("--csv", action="store", type=str, dest="in_csv")
    parser.add_argument("--reads", action="store", type=str, dest="reads")

    parser.add_argument("--bin_factor", action="store", type=str, dest="bin_factor")
    parser.add_argument("--references", action="store", type=str, dest="references")

    parser.add_argument("--sample", action="store", type=str, dest="sample")

    parser.add_argument("--min_reads", action="store", type=int, dest="min_reads")
    parser.add_argument("--min_pcent", action="store", type=float, dest="min_pcent")

    parser.add_argument("--output_path", action="store", type=str, dest="output_path")

    return parser.parse_args()


def parse_csv(file):
    reader = csv.DictReader(file)
    data = [r for r in reader]
    return data

def get_files_and_stems(in_csv, bin_factor, min_reads, output_path, references, sample):

    bin_factor_counts = {
        "Amp123": collections.Counter(),
        "Amp4": collections.Counter(),
        "Amp5": collections.Counter()
        }
        
    best_ref_counts = {
            "Amp123": collections.Counter(),
            "Amp4": collections.Counter(),
            "Amp5": collections.Counter()
        }

    bin_factor_map= {"Amp123": {},"Amp4": {},"Amp5": {}}
    min_reads_per_amp = {"Amp123":min_reads*3,"Amp4":min_reads,"Amp5":min_reads}

    with open(in_csv,"r") as f:
        for row in parse_csv(f):
            if row["best_reference"] in ['*','?','none']:
                pass
            else:
                bin_factor_counts[row["amplicon"]][row[bin_factor]]+=1 #e.g. genotype count per amplicon
                bin_factor_map[row["amplicon"]][row["best_reference"]]=row[bin_factor] #e.g. map from accession no. to genotype
                best_ref_counts[row["amplicon"]][row["best_reference"]]+=1 #e.g. accession no. count per amplicon

    greater_than_minimum_bin_factors = defaultdict(list)
    relevant_accessions = defaultdict(list)

    for amplicon in ["Amp123","Amp4","Amp5"]:
        
        for factor in bin_factor_counts[amplicon]:
            if bin_factor_counts[amplicon][factor] > min_reads_per_amp[amplicon]: #amp123 requires 3 times as many reads
                greater_than_minimum_bin_factors[amplicon].append(factor)
        
        for factor in greater_than_minimum_bin_factors[amplicon]: #GII16
            for best_ref in best_ref_counts[amplicon]:#accession nos
                
                if factor == bin_factor_map[amplicon][best_ref]: #if GII16 == the loc_genotype for that accession no
                    relevant_accessions[amplicon].append(best_ref)
    
    common_to_all = list(set(relevant_accessions["Amp123"]) & set(relevant_accessions["Amp4"]) & set(relevant_accessions["Amp5"]))

    analysis_guides = defaultdict(list)

    if common_to_all:
        for accession in common_to_all:
            for amplicon in ["Amp123","Amp4","Amp5"]:
                info_to_pass = (bin_factor_map[amplicon][accession],accession)
                analysis_guides[amplicon].append(info_to_pass)
                
    else:
        for amplicon in ["Amp123","Amp4","Amp5"]:
            for factor in greater_than_minimum_bin_factors[amplicon]:
                    
                    top_ref_per_factor = []
                    for i in best_ref_counts[amplicon]:
                        if bin_factor_map[amplicon][i] == factor:
                            top_ref_per_factor.append((i,best_ref_counts[amplicon][i]))
                    top_ref_per_factor = sorted(top_ref_per_factor, key= lambda x : x[1], reverse=True)[0]

                    info_to_pass = (factor,top_ref_per_factor[0])
                    analysis_guides[amplicon].append(info_to_pass)
    
    read_dict = defaultdict(list)
    index = 0
    with open(in_csv,"r") as f:
        for row in parse_csv(f):
            index+=1
            if row["best_reference"] in ['*','?','none']:
                pass
            else:
                amp = row["amplicon"] 
                factor = row[bin_factor]
                orf = ""
                if amp in ["Amp123","Amp5"]:
                    if amp == "Amp123":
                        orf = "ORF1"
                    else:
                        orf = "ORF23"
                    for guide in analysis_guides[amp]:
                        stem = f"{orf}_{guide[0]}"
                        if factor==guide[0]:
                            read_dict[stem].append(index)

                        write_reference_file(output_path, stem, sample, guide[1], references)
                        
                elif amp == "Amp4":
                    for guide in analysis_guides["Amp123"]:
                        stem1 = f"ORF1_{guide[0]}"
                        read_dict[stem1].append(index)
                    for guide in analysis_guides["Amp5"]:
                        stem2 = f"ORF23_{guide[0]}"
                        read_dict[stem2].append(index)


    return read_dict

def file_writer(read_dict,reads,output_path,in_csv,sample):
    analysis_stems = ''
    for stem in read_dict:
        analysis_stem = f"{sample}_{stem}"
        analysis_stems+=analysis_stem+','
        new_seq_file = output_path + f"/{analysis_stem}.fastq"
        new_csv_file = output_path + f"/{analysis_stem}.csv"
        record_index = 0
        csv_index = -1
        records = []
        with open(new_seq_file,"w") as fw:
            for record in SeqIO.parse(reads,"fastq"):
                record_index +=1
                if record_index in read_dict[stem]:
                    records.append(record)
            SeqIO.write(records, fw, "fastq")
        with open(new_csv_file,"w") as fw:
            with open(in_csv,"r") as f:
                for l in f:
                    l = l.rstrip('\n')
                    csv_index +=1
                    if csv_index ==0:
                        fw.write(l + '\n')
                    if csv_index in read_dict[stem]:
                        fw.write(l + '\n')
    print(analysis_stems.rstrip(','))

def write_reference_file(output_path, stem, sample, accession, references):
    filename = f'{output_path}/{sample}_{stem}.fasta'
    with open(filename,"w") as fw:
        for record in SeqIO.parse(references,"fasta"):
            if record.id == accession:
                fw.write(f">{record.description}\n{record.seq}\n")

if __name__ == '__main__':

    args = parse_args()

    read_dict = get_files_and_stems(args.in_csv, args.bin_factor, args.min_reads, args.output_path, args.references, args.sample)

    file_writer(read_dict,args.reads,args.output_path,args.in_csv,args.sample)
