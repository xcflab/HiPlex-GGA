from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import json
from io import StringIO
import pandas as pd
import os

# 定义四个恒定区
CONSTANT_REGION_1 = "aagcttatgagcTCTAGA".upper()  # 
CONSTANT_REGION_2 = "ATCCCGTCGAGTTTGACCACGTG".upper()  # 
CONSTANT_REGION_3 = "cgggattcgaagaCCATGG".upper()  # 
CONSTANT_REGION_4 = "CATATG".upper()  # NdeI
pattern_1 = CONSTANT_REGION_1.upper() + "." * 19 + CONSTANT_REGION_2 
# 条形码的预期长度
BARCODE_LENGTH = 19

# 批次的长度和三种允许的序列
BATCH_LENGTH = 6
ALLOWED_BATCHES = ["atgagc", "tcaacg", "ggaaac"]  # 

# 定义输入和输出文件路径
#input_fasta_file = "/xcfhome/zewu/ngs/N2422148_80-1626558504_2024-09-11/240909-A00599B/test.fasta"
input_fasta_file = "/xcfhome/zewu/ngs/N2422148_80-1626558504_2024-09-11/240909-A00599B/PE150-ZSZ-0903-LGJ3824_trimmed_merged.fasta"

def reverse_fasta(record):
    """
    read and reverse sequences in fasta file
    """
    reversed_seq = record.seq.reverse_complement()  # reverse the sequence
    reversed_record = record  # copy the record
    reversed_record.seq = reversed_seq  # replace the sequence with reversed sequence
    return reversed_record

def find_constant_regions(seq):
    """找到恒定区，并检查是否符合要求"""
    try:
        # 查找四个恒定区的位置
        start1 = seq.index(CONSTANT_REGION_1)
        start2 = seq.index(CONSTANT_REGION_2, start1 + len(CONSTANT_REGION_1))
        start3 = seq.index(CONSTANT_REGION_3, start2 + len(CONSTANT_REGION_2))
        start4 = seq.index(CONSTANT_REGION_4, start3 + len(CONSTANT_REGION_3))
        
        # 提取条形码
        barcode = seq[start1 + len(CONSTANT_REGION_1):start2]
        if len(barcode) != BARCODE_LENGTH:
            return None, None
        
        # 提取批次
        batch = seq[start3 + len(CONSTANT_REGION_3):start4]
        if len(batch) != BATCH_LENGTH or batch not in ALLOWED_BATCHES:
            return None, None
        
        return barcode, batch
    except ValueError:
        # 如果找不到任何恒定区，返回None
        return None, None
    
def process_fasta(input_file):
    """
    Find the sequences that contain subpool sequence and filter the sequences with length 227
    """
    batch_1_records = []
    pattern_2_1 = CONSTANT_REGION_3.upper() + ALLOWED_BATCHES[0].upper() + CONSTANT_REGION_4.upper()
    batch_2_records = []
    pattern_2_2 = CONSTANT_REGION_3.upper() + ALLOWED_BATCHES[1].upper() + CONSTANT_REGION_4.upper()
    batch_3_records = []
    pattern_2_3 = CONSTANT_REGION_3.upper() + ALLOWED_BATCHES[2].upper() + CONSTANT_REGION_4.upper()
    
    for record in SeqIO.parse(input_file, "fasta"):
        if re.search(pattern_1, str(record.seq)) and len(record.seq) == 227:
            if re.search(pattern_2_1, str(record.seq)):
                batch_1_records.append(record)
            elif re.search(pattern_2_2, str(record.seq)):
                batch_2_records.append(record)
            elif re.search(pattern_2_3, str(record.seq)):
                batch_3_records.append(record)
        else:
            record = reverse_fasta(record)
            if re.search(pattern_1, str(record.seq)) and len(record.seq) == 227:
                if re.search(pattern_2_1, str(record.seq)):
                    batch_1_records.append(record)
                elif re.search(pattern_2_2, str(record.seq)):
                    batch_2_records.append(record)
                elif re.search(pattern_2_3, str(record.seq)):
                    batch_3_records.append(record)
                  
    batch = [batch_1_records, batch_2_records, batch_3_records]   
    return batch

def count_reads(records):
    """
    
    """
    total_records = {}
    for record in records:
        seq = record.seq
        umi_start = seq.find(CONSTANT_REGION_1) + len(CONSTANT_REGION_1)
        umi_end = seq.find(CONSTANT_REGION_2)
        umi = seq[umi_start:umi_end]
        total_records[umi] = total_records.get(umi, 0) + 1 # count the number of reads with the same umi
    sorted_records = sorted(total_records.items(), key=lambda x: x[1], reverse=True)
    return sorted_records
    
def sta_data(input_file):
    """find the sequences that contain subpool sequence and filter the sequences with length larger than 1500"""
    
    batch = process_fasta(input_file)
    for i in range(3):
        print(f"批次 {i+1} 共有 {len(batch[i])} 条序列。")
        count = count_reads(batch[i])
        df = pd.DataFrame(count, columns=["UMI", "Reads"])
        csv_file = f"/xcfhome/zewu/ngs/N2422148_80-1626558504_2024-09-11/240909-A00599B/count_batch{i+1}.csv"
        df.to_csv(csv_file, index=False)
    
if __name__ == "__main__":
    sta_data(input_fasta_file)
    print("序列处理完成。")