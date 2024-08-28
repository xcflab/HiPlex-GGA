from pyfasta import Fasta
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import json
from io import StringIO
import pandas as pd
import os

# read fasta file
work_path = "/xcfhome/zewu/pacbio" #change to your own path
test_file = "test.fasta"
official_file = "zsz_pacbio_240605l.fasta"
pacbio_fasta_file = os.path.join(work_path, official_file)  # fasta file path zsz_pacbio_240605l.fasta
template_fasta_file = "sub31_complete.fasta"
template_fasta_file_path = os.path.join(work_path, template_fasta_file)  # template fasta file path
split_template_fasta_file = "sub31_split.fasta"
split_template_fasta_file_path = os.path.join(work_path, split_template_fasta_file) # split template fasta file path
sequence_barcode_file = "sequence_barcode_list.txt"
sequence_barcode_file_path = os.path.join(work_path, sequence_barcode_file) # output fasta file path
batch_1 = "AAGCTTatgagcTCTAGANNNNNNNNNNNNNNNNNNNATCccgtcgagtttgaCCACGTG"  # batch 1 sequence
batch_2 = "AAGCTTtcaacgNNNNNNNNNNNNNNNCCATGG"  # batch 2 sequence

def reverse_fasta(record):
    # read and reverse every sequences in fasta file
    reversed_seq = record.seq.reverse_complement()  # reverse the sequence
    reversed_record = record  # copy the record
    reversed_record.seq = reversed_seq  # replace the sequence with reversed sequence
    return reversed_record

def adjust_record(record, start, end):
    # adjust the record with start and end position
    agjusted_seq = record.seq[start:end]
    adjusted_record = record
    adjusted_record.seq = Seq(agjusted_seq)
    return adjusted_record

def grep_subpool_sequence(fasta, subpool_sequence, length=900):
    # find the sequences that contain subpool sequence and filter the sequences with length larger than 1500
    subpool_records = []  
    pattern = subpool_sequence.replace("N", ".")
    pattern = pattern.upper()
    for record in fasta:
        if re.search(pattern, str(record.seq)) and len(record.seq) < length:
            subpool_records.append(record)  # append the record to the list
        else:
            record = reverse_fasta(record)
            if re.search(pattern, str(record.seq)) and len(record.seq) < length:
                subpool_records.append(record)
    return subpool_records

def sort_umi(records, cons_1, cons_2, complete_sequences_set):
    HindIII = "AAGCTT" # HindIII site
    EcoRI = "GAATTC" # EcoRI site
    umis = {} # store the umi, unique molecular identifier
    
    for record in records:
        umi = record.seq[record.seq.find(cons_1) + len(cons_1) : record.seq.find(cons_2)]
        begin = record.seq.find(HindIII) + len(HindIII)
        end = record.seq.find(EcoRI)
        sequence = record.seq[begin:end]
        
        if umi in umis.keys():
            record = SeqRecord(sequence, record.id, description=str(umi))
            umis[umi].append(record)
        else:
            record = SeqRecord(sequence, record.id, description=str(umi))
            umis[umi] = [record]
            
    return umis

def check_unique_reads(reads_list):
    seen = set()
    for reads in reads_list:
        if reads.seq not in seen:
            seen.add(reads.seq)
            
    if len(seen) > 2: # change here for the number of unique reads
        return True
    else:
        return False

def umi_statistics(umi, reads_list, complete_sequences_set):
    # 
    ef_reads = []
    for read in reads_list:
        if any(seq in read.seq for seq in complete_sequences_set):
            ef_reads.append(read)
    
    if check_unique_reads(reads_list):
        umi_file = str(umi) + ".fasta"
        umi_file_path = os.path.join(work_path, "umi", umi_file)
        with open(umi_file_path, "w") as output_handle:
            SeqIO.write(reads_list, output_handle, "fasta")
    else:
        umi_file_path = None
    
    num_reads = len(reads_list)
    num_ef_reads = len(ef_reads)
    umi_dict={"total_reads": num_reads, "ef_reads": num_ef_reads, "umi_file_path": umi_file_path}
    
    return umi_dict

def count_all_sequences(records, complete_sequences):
    cons_1 = "ATGAGCTCTAGA" # constant region 1 -> batch1+XbaI
    cons_2 = "ATCCCGTCGAGTTTGACCACGTG" # constant region 2 -> AdapterF
    
    complete_sequences_set = set(sequence.seq for sequence in complete_sequences)
    umis = sort_umi(records, cons_1, cons_2, complete_sequences_set)
    
    total_error_free = 0
    umi_sta = {}
    for key,value in umis.items():
        umi_sta[key] = umi_statistics(key, value, complete_sequences_set)
        if umi_sta[key]["ef_reads"] > 0:
            total_error_free += 1
    
    total_umi = len(umis)
    ratio = total_error_free / total_umi
    # store the umi statistics to csv file
    df = pd.DataFrame.from_dict(umi_sta, orient='index')
    csv_file = os.path.join(work_path, "umi_statistics.csv")
    df.to_csv(csv_file)
    
    return [total_umi, total_error_free, ratio]

##### MAIN #####
# find the sequences that contain batch 1 sequence
pacbio_fasta = SeqIO.parse(pacbio_fasta_file, "fasta")
real_subpool_records_1 = grep_subpool_sequence(pacbio_fasta, batch_1)
# read template fasta file
template_sequences = list(SeqIO.parse(template_fasta_file_path, "fasta"))
# count the number of sequences that contain the template sequence (batch 1)
error_free_batch1 = count_all_sequences(real_subpool_records_1, template_sequences)
print("total reads: {}, umi: {}, error-free umi: {}, error-free ratio: {:.2f}".format(len(real_subpool_records_1), *error_free_batch1))