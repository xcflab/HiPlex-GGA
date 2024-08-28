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
pacbio_fasta_file = os.path.join(work_path, test_file)  # fasta file path zsz_pacbio_240605l.fasta
template_fasta_file = "/xcfhome/zewu/pacbio/sub31_complete.fasta"  # template fasta file path
split_template_fasta_file = "/xcfhome/zewu/pacbio/sub31_split.fasta"  # split template fasta file path
#output_fasta = "/xcfhome/zewu/pacbio/sub31_seq1_filtered.fasta"  # output fasta file path
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
    for record in SeqIO.parse(pacbio_fasta_file, "fasta"):
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
    if len(seen) > 1:
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

def count_single_sequence(records, sequence, sequence_barcode, umi_length=15, bamhi_ecori="GGATCCgaattc"):
    sequence_barcode = sequence_barcode.upper()
    bamhi_ecori = bamhi_ecori.upper()
    # count the number of sequences that contain the sequence
    total_umis = set() # store the umi, unique molecular identifier
    ef_umis = set()
    total_matching_records = []
    ef_matching_records = []
    for record in records:
        if sequence_barcode in record.seq:
            barcode_pos = record.seq.find(sequence_barcode)
            bamhi_ecori_pos = record.seq.find(bamhi_ecori) + len(bamhi_ecori)
            if barcode_pos != -1:
                umi_start = barcode_pos - umi_length
                umi_end = barcode_pos
                umi = record.seq[umi_start:umi_end]
                if umi not in total_umis:
                    total_umis.add(umi)
                    adjusted_record = adjust_record(record, umi_start - 12, bamhi_ecori_pos)
                    total_matching_records.append(adjusted_record)
                if sequence in record.seq and umi not in ef_umis:
                    ef_umis.add(umi)
                    ef_matching_records.append(adjusted_record)

    with open("/xcfhome/zewu/pacbio/seq4_total_umi.fasta", "w") as output_handle:
        SeqIO.write(total_matching_records, output_handle, "fasta")
    with open("/xcfhome/zewu/pacbio/seq4_ef_umi.fasta", "w") as output_handle:
        SeqIO.write(ef_matching_records, output_handle, "fasta")
    total_error_free = len(ef_umis)
    count_total = len(total_umis)
    ratio = total_error_free / count_total
    
    return [count_total, total_error_free, ratio]

def count_all_sequences(records, complete_sequences, split_sequence):
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

def store_all_sequences(records, complete_sequences, split_sequence):
    cons_1 = "ATGAGCTCTAGA" # constant region 1 -> batch1+XbaI
    cons_2 = "ATCCCGTCGAGTTTGACCACGTG" # constant region 2 -> AdapterF
    complete_sequences_dict = {sequence.seq: sequence for sequence in complete_sequences}
    # count the number of sequences that contain the sequence
    umis = {} # store the umi, unique molecular identifier
    for record in records:
        umi = record.seq[record.seq.find(cons_1) + len(cons_1) : record.seq.find(cons_2)]
        begin = record.seq.find(HindIII) + len(HindIII)
        end = record.seq.find(EcoRI)
        sequence = record.seq[begin:end]
        if umi in umis.keys():
            record = [sequence, record.id]
            umis[umi].append(record)
        else:
            record = [sequence, record.id]
            umis[umi] = [record]
    
    check_umis = []
    for key, value in umis.items():
        if len(value) > 7:
            description = str(key)
            output = StringIO()
            for record in value:
                seq_record = SeqRecord(record[0], record[1], description=description)
                check_umis.append(seq_record)
                
    output_fasta = "/xcfhome/zewu/pacbio/more_than_2_umi.fasta"
    with open(output_fasta, "w") as output_handle:
        SeqIO.write(check_umis, output_handle, "fasta")

# pacbio_fasta = reverse_fasta(pacbio_fasta_file)  # read pacbio fasta file
# find the sequences that contain batch 1 sequence
real_subpool_records_1 = grep_subpool_sequence(pacbio_fasta_file, batch_1)
real_subpool_records_2 = grep_subpool_sequence(pacbio_fasta_file, batch_2)
# read template fasta file
template_sequences = list(SeqIO.parse(template_fasta_file, "fasta"))
frag_sequences = list(SeqIO.parse(split_template_fasta_file, "fasta"))

# count the number of sequences that contain the template sequence (batch 1)
error_free_batch1 = count_all_sequences(real_subpool_records_1, template_sequences, frag_sequences)
# print("total reads: {}, umi: {}, error-free umi: {}, error-free ratio: {:.2f}".format(len(real_subpool_records_1), *error_free_batch1))
# store_umi = store_all_sequences(real_subpool_records_1, template_sequences, frag_sequences)
# count the number of sequences that contain the template sequence
# temp = template_sequences[3]
# sequence_barcode = "CCATGG" + "ggcgactccataa" + "CATATG"
# test = count_single_sequence(real_subpool_records_2, temp.seq, sequence_barcode)
# print(*test)
