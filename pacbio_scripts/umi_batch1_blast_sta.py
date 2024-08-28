#from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import json
from io import StringIO
import pandas as pd
import os
import subprocess
    
# read fasta file
work_path = "/xcfhome/zewu/pacbio" #change to your own path
test_file = "test.fasta"
official_file = "zsz_pacbio_240605l.fasta"
pacbio_fasta_file = os.path.join(work_path, test_file)  # fasta file path zsz_pacbio_240605l.fasta
template_fasta_file = "sub31_complete.fasta"
template_fasta_file_path = os.path.join(work_path, template_fasta_file)  # template fasta file path
split_template_fasta_file = "sub31_split.fasta"
split_template_fasta_file_path = os.path.join(work_path, split_template_fasta_file) # split template fasta file path
sequence_barcode_file = "sequence_barcode_list.txt"
sequence_barcode_file_path = os.path.join(work_path, sequence_barcode_file) # output fasta file path
batch_1 = "AAGCTTatgagcTCTAGANNNNNNNNNNNNNNNNNNNATCccgtcgagtttgaCCACGTG"  # batch 1 sequence
batch_2 = "AAGCTTtcaacgNNNNNNNNNNNNNNNCCATGG"  # batch 2 sequence
subpool_sequence_begin = "AAGCTTatgagcTCTAGANNNNNNNNNNNNNNNNNNNATCccgtcgagtttgaCCACGTG"
subpool_sequence_end = "AAGCTTtcaacgNNNNNNNNNNNNNNNCCATGG"

#enzyme_list
HindIII = 'AAGCTT'
NcoI = 'CCATGG'
NdeI = 'CATATG'
XhoI = 'CTCGAG'
BamHI = 'GGATCC'
EcoRI = 'GAATTC'
#input the site of enzyme
#5'
BsaI = 'GGTCTCA'
BsmBI = 'CGTCTCA'
#3'
bsai = 'AGAGACC'
bsmbi = 'AGAGACG'

def read_sequence_barcode(infile, seq_num):
    """
    Read the sequence barcode file, obetaining the correct sequences from the expected format.
    """
    handle = open(infile, 'r')
    # get all lines, skip final newline and drop first line that has titles:
    lines = [ line[:-1] for line in handle.readlines() ][1:]
    seq_barcode_list=[]
    for i in range(seq_num):
        l = lines[i].split()[1:]
        seq_barcode_list.append(l)
    return seq_barcode_list

def umi_correction(umi, reads_list):
    corrected = umi_tools.collapse(reads_list)
    return corrected

def find_sequence_belong_gene(seq_name, seq):
    """
    Find the sequence that belongs to the gene
    """
    seq = seq.upper()
    return seq

'''
def blast_bio(seq, db):
    # 配置命令行参数
    blastn_cline = NcbiblastnCommandline(
        query=seq, 
        db=template_fasta_file_path, 
        evalue=0.001, 
        outfmt=6, 
        out="output.txt"
    )

    # 运行 BLAST
    stdout, stderr = blastn_cline()

    # 读取并打印输出文件
    with open("output.txt") as result_file:
        print(result_file.read())
'''

def blast_command(umi, seq, db):
    # create seq file
    file_name = str(umi) + ".fasta"
    seq_path = os.path.join(work_path, file_name)
    with open(seq_path, "w") as seq_file:
        SeqIO.write(seq, seq_file, "fasta")
    
    temp_file = os.path.join(work_path, "output.txt.tmp")
    # define the blastn command
    blastn_cmd = [
        "blastn", 
        "-query", seq_path, 
        "-db", template_fasta_file_path, 
        "-out", temp_file, 
        "-evalue", "1e-5", 
        "-outfmt", "6"
    ]
    
    # run the blastn command
    subprocess.run(blastn_cmd)

    # read and print the result
    with open(temp_file, "r") as result_file:
        blast_list = []
        for l in result_file.readlines():
            blast_result = l.split()
        
        print(blast_result)
        
    #delete the temp file
    os.remove(seq_path)
    os.remove(temp_file)

def assembly(infile, frag_num=4):
    """
    Read the split template fasta file and get the sequences
    """
    handle = open(infile, 'r')
    # get all lines, skip final newline and drop first line that has titles:
    split_record_list = list(line for line in handle.readlines())
    seq_dict = {}
    for i in range(0, len(split_record_list), 2 * frag_num):
        for j in range(0, frag_num*2, 2):
            key = split_record_list[i+j]
            seq = split_record_list[i+j+1]
            gene_seq = find_sequence_belong_gene(key, seq)
            seq_dict[key] = gene_seq
            
    return seq_dict

def reverse_fasta(record):
    """
    read and reverse sequences
    """
    reversed_seq = record.seq.reverse_complement()  # reverse the sequence
    reversed_record = record  # copy the record
    reversed_record.seq = reversed_seq  # replace the sequence with reversed sequence
    return reversed_record

def adjust_record(record, start, end):
    """
    Adjust the record with start and end position
    """
    agjusted_seq = record.seq[start:end]
    adjusted_record = record
    adjusted_record.seq = Seq(agjusted_seq)
    return adjusted_record

def grep_subpool_sequence(fasta, subpool_sequence, length=900):
    """
    Find the sequences that contain subpool sequence and keep the sequences with length shorter than 900
    """
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
    """
    Read every umi and statistic the reads
    """
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
    
    blast_command(umi, reads_list, template_fasta_file_path)
    
    num_reads = len(reads_list)
    num_ef_reads = len(ef_reads)
    umi_dict={"total_reads": num_reads, "ef_reads": num_ef_reads, "umi_file_path": umi_file_path}
    
    return umi_dict

def seq_statistics(umi, reads_list, complete_sequences_set):
    """
    Read every umi and statistic the reads
    """
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
    
    blast_command(umi, reads_list, template_fasta_file_path)
    
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
        # umi_sta[key] = umi_statistics(key, value, complete_sequences_set)
        # umi_sta[key] = seq_statistics(key, value, complete_sequences_set)
        umi_sta[key] = umi_correction(key, value)
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
# blast_command(pacbio_fasta_file, template_fasta_file_path)
real_subpool_records_1 = grep_subpool_sequence(pacbio_fasta, batch_1)
# read template fasta file
template_sequences = list(SeqIO.parse(template_fasta_file_path, "fasta"))
seq_barcode = read_sequence_barcode(sequence_barcode_file_path, 32)
split_template_sequences = assembly(split_template_fasta_file_path)
# count the number of sequences that contain the template sequence (batch 1)
error_free_batch1 = count_all_sequences(real_subpool_records_1, template_sequences)
print("total reads: {}, umi: {}, error-free umi: {}, error-free ratio: {:.2f}".format(len(real_subpool_records_1), *error_free_batch1))