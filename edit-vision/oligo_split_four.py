#!/usr/bin/env python
# Mar 12, 2024
# Zhien Wu modified
import numpy as np
from tqdm import tqdm
import sys
import random
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

#overhand_list
#The best overhangs are included in overhangs_1st. The second good overhangs are included in overhangs_2nd.
overhangs_1st = ['TCAA', 'ATAA', 'TAGA', 'GTTA', 'AGGG', 'CCTA', 'AAGA', 'TCCA', 'AGAA', 'AAAG', 'ACAT', 'GGGA', 'GCAA', 'TAAA', 'TGAA', 'GACA', 'ACGA', 'GGTA', 'ATCC', 'ATTG', 'CTAC', 'AAGC', 'CATC', 'ACTC', 'CACA', 'CTAA', 'GAAA', 'AGAC', 'AGCA', 'CGTA', 'ACCG', 'CAGA', 'AACT', 'AATA', 'GCAC', 'CCAG', 'CAAG', 'AAAT', 'ATCA', 'CAGG', 'CATA', 'GGAA', 'AGGA', 'ACGC', 'ATAC', 'CTCA', 'GCCA', 'CCGA', 'ACAG', 'AATC', 'CAGC', 'AAAA', 'AGTG', 'CGCA', 'AACG', 'GAGA', 'ACTA', 'TACA', 'ATGA', 'CGAC', 'CGAA', 'AGCC']
overhangs_2nd = ['GGGG', 'CCCC', 'GGGC', 'GCCC', 'GGCG', 'CGCC', 'GCGG', 'CCGC', 'CGGG', 'CCCG', 'GGCC', 'CCGG', 'GCCG', 'CGGC', 'CGCG', 'GCGC']#too much?

#enzyme_list
#HindIII = 'AAGCTT'
NcoI = 'CCATGG'
NdeI = 'CATATG'
XhoI = 'CTCGAG'
BamHI = 'GGATCC'
#EcoRI = 'GAATTC'
#input the site of enzyme
#3'
bsai = 'AGAGACC'
bsmi = 'AGAGACG'
#5'
BsaI = 'GGTCTCA'
BsmI = 'CGTCTCA'

#function:
def shuffle_dict(input_dict):
    keys = list(input_dict.keys())
    random.shuffle(keys)
    shuffled_dict = {key: input_dict[key] for key in keys}
    return shuffled_dict

def read_adapters(infile, frag_num):
    """
    Reads in the adapter file, obetaining the correct sequences from the expected format.
    """
    handle = open(infile, 'r')
    # get all lines, skip final newline and drop first line that has titles:
    lines = [ line[:-1] for line in handle.readlines() ][1:]
    adapter_list=[]
    for i in range(frag_num * 2):
        l = [ line.split()[i+1] for line in lines ]
        adapter_list.append(l)
    return adapter_list

def read_sequence_barcode(infile, seq_num):
    """
    Reads in the sequence barcode file, obetaining the correct sequences from the expected format.
    """
    handle = open(infile, 'r')
    # get all lines, skip final newline and drop first line that has titles:
    lines = [ line[:-1] for line in handle.readlines() ][1:]
    seq_barcode_list=[]
    for i in range(seq_num):
        l = lines[i].split()[1:]
        seq_barcode_list.append(l)
    return seq_barcode_list

def count_frag_lenth(seq_len, frag_th, frag_num_left, max_inner_len=220, max_last_len=210, max_fist_len=210, min_len=110):
    frag_min_len = seq_len - (max_inner_len - 4) * ( frag_num_left - 1 ) - ( max_last_len - 4)
    if frag_min_len < min_len:
        frag_min_len = min_len
    
    if frag_th == 1:
        assert frag_min_len <= max_fist_len
    else:
        assert frag_min_len <= max_inner_len

    if frag_num_left == 1:
        last_len = seq_len - max_inner_len
        if last_len < min_len:
            max_inner_len = seq_len - min_len * ( frag_num_left - 1 )
        
    return frag_min_len

def check_overhang_unique(end_list, to_check_end):
    # 检查粘性末端是否唯一
    flag = True
    if to_check_end in end_list:
        flag = False
    return flag

def split_with_pos(sequence, pos, overhang_len=4):
    end_len_1 = overhang_len
    seq_1 = sequence[0:pos]
    seq_2 = sequence[pos-end_len_1:]
    overhang = Seq(sequence[pos-end_len_1:pos])
    overhang_pair = overhang.reverse_complement()
    return seq_1, seq_2, overhang, overhang_pair

def split_with_pos_re(sequence, pos, overhang_len=4):
    end_len_1 = overhang_len
    seq_1 = sequence[pos:]
    seq_2 = sequence[0:pos-end_len_1]
    overhang = Seq(sequence[pos-end_len_1:pos])
    overhang_pair = overhang.reverse_complement()
    return seq_1, seq_2, overhang, overhang_pair

def find_overhang_from_list_max(end_pos, i, max_len_5, max_inner_len, overhang_list, seq, overhang_len):
    (frag, rest_seq, overhang, overhang_pair) =  (None, None, None, None)
    split_successful_tag = False
    
    max_len = max_len_5 if i == 1 else max_inner_len
    
    for j in range(end_pos, max_len + 1):
        frag, rest_seq, overhang, overhang_pair = split_with_pos(seq, j, overhang_len)
        if overhang in overhangs_1st and check_overhang_unique(overhang_list, overhang):
            split_successful_tag = True
            break
        
    if split_successful_tag:
        return frag, rest_seq, overhang, overhang_pair
    else:
        return split_successful_tag

def find_overhang_from_list_max_re(end_pos, i, max_len_3, max_inner_len, overhang_list, seq, overhang_len, frag_num):
    (frag, rest_seq, overhang, overhang_pair) =  (None, None, None, None)
    split_successful_tag = False
    
    max_len = len(seq) - max_len_3 if i == 1 else len(seq) - max_inner_len
 
    for j in range(end_pos + 1, max_len, -1 ):
        frag, rest_seq, overhang, overhang_pair = split_with_pos_re(seq, j, overhang_len)
        if overhang in overhangs_1st and check_overhang_unique(overhang_list, overhang):
            split_successful_tag = True
            break
        
    if split_successful_tag:
        return frag, rest_seq, overhang, overhang_pair
    else:
        return split_successful_tag
    
def find_overhang_from_list_min(end_pos, frag_min_len, overhang_list, seq, overhang_len):
    (frag, rest_seq, overhang, overhang_pair) =  (None, None, None, None)
    split_successful_tag = False
    
    for j in range(end_pos, frag_min_len - 1, -1):
        frag, rest_seq, overhang, overhang_pair = split_with_pos(seq, j, overhang_len)
        if overhang in overhangs_1st and check_overhang_unique(overhang_list, overhang):
            split_successful_tag = True
            break
        
    if split_successful_tag:
        return frag, rest_seq, overhang, overhang_pair
    else:
        return split_successful_tag

def find_overhang_from_list_min_re(end_pos, frag_min_len, overhang_list, seq, overhang_len):
    (frag, rest_seq, overhang, overhang_pair) =  (None, None, None, None)
    split_successful_tag = False
    
    min_len = len(seq) - frag_min_len
    for j in range(end_pos, min_len + 1):
        frag, rest_seq, overhang, overhang_pair = split_with_pos_re(seq, j, overhang_len)
        if overhang in overhangs_1st and check_overhang_unique(overhang_list, overhang):
            split_successful_tag = True
            break
        
    if split_successful_tag:
        return frag, rest_seq, overhang, overhang_pair
    else:
        return split_successful_tag
def find_overhang_from_2nd(end_pos, i, max_len_5, max_inner_len, frag_min_len, overhang_list, seq, overhang_len):
    (frag, rest_seq, overhang, overhang_pair) =  (None, None, None, None)
    split_successful_tag = False
    
    max_len = max_len_5 if i == 1 else max_inner_len
    
    for j in range(end_pos, max_len + 1):
        frag, rest_seq, overhang, overhang_pair = split_with_pos(seq, j, overhang_len)
        if overhang not in overhangs_2nd and check_overhang_unique(overhang_list, overhang):
            split_successful_tag = True
            break
    while overhang in overhangs_2nd or not check_overhang_unique(overhang_list, overhang):
        for j in range(frag_min_len, end_pos):
            frag, rest_seq, overhang, overhang_pair = split_with_pos(seq, j, overhang_len)
            if overhang not in overhangs_2nd and check_overhang_unique(overhang_list, overhang):
                split_successful_tag = True
                break
            
    if split_successful_tag:
        return frag, rest_seq, overhang, overhang_pair
    else:
        return split_successful_tag

def find_overhang_from_2nd_re(end_pos, i, max_len_3, max_inner_len, frag_min_len, overhang_list, seq, overhang_len, frag_num):
    (frag, rest_seq, overhang, overhang_pair) =  (None, None, None, None)
    split_successful_tag = False
    
    max_len = len(seq) - max_len_3 if i == 1 else len(seq) - max_inner_len
    
    for j in range(end_pos + 1, max_len, -1):
            frag, rest_seq, overhang, overhang_pair = split_with_pos_re(seq, j, overhang_len)
            if overhang not in overhangs_2nd and check_overhang_unique(overhang_list, overhang):
                split_successful_tag = True
                break
    min_len = len(seq) - frag_min_len
    while overhang in overhangs_2nd or not check_overhang_unique(overhang_list, overhang):
        for j in range(end_pos, min_len + 1):
            frag, rest_seq, overhang, overhang_pair = split_with_pos_re(seq, j, overhang_len)
            if overhang not in overhangs_2nd and check_overhang_unique(overhang_list, overhang):
                split_successful_tag = True
                break
            
    if split_successful_tag:
        return frag, rest_seq, overhang, overhang_pair
    else:
        return split_successful_tag

def find_unique_overhang(end_pos, i, max_len_5, max_len_3, max_inner_len, frag_min_len, overhang_list, seq, seq_num, frag_num, overhang_len = 4):
    split_successful_tag = False
    while not split_successful_tag:
        if seq_num % 2 == 0:
            #Choose overhang from end_pos to max_len
            result = find_overhang_from_list_max(end_pos, i, max_len_5, max_inner_len, overhang_list, seq, overhang_len)
            try_times = 0
            while isinstance(result, bool):
                result = find_unique_overhang(end_pos, i, max_len_5, max_len_3, max_inner_len, frag_min_len, overhang_list[i-1], seq, seq_num, frag_num)
                try_times += 1
                if try_times == 10:
                    print("try too many times")
                    break
            else:
                frag, rest_seq, overhang, overhang_pair = result
            #If there is no unique overhang from end_pos to max_len, choose from min_len to end_pos
            if overhang not in overhangs_1st or not check_overhang_unique(overhang_list, overhang):
                result = find_overhang_from_list_min(end_pos, frag_min_len, overhang_list, seq, overhang_len)
                if isinstance(result, bool):
                    # print("Can't split sequence, try again!")
                    split_successful_tag = False
                    break
                else:
                    frag, rest_seq, overhang, overhang_pair = result
                    split_successful_tag = True

            if overhang not in overhangs_1st or not check_overhang_unique(overhang_list, overhang):
                frag, rest_seq, overhang, overhang_pair = find_overhang_from_2nd(end_pos, i, max_len_5, max_inner_len, frag_min_len, overhang_list, seq, overhang_len)
            
        #  if overhang in overhangs_1st:
        #      print("choose from best group: ", overhang)
        #  else:
        #      print("overhang:", overhang )
        #  if not check_overhang_unique(overhang_list, overhang):
        #      print("no unique")
                
        else:
            #Choose overhang from end_pos to max_len
            frag, rest_seq, overhang, overhang_pair = find_overhang_from_list_max_re(end_pos, i, max_len_3, max_inner_len, overhang_list, seq, overhang_len, frag_num)

            #If there is no unique overhang from end_pos to max_len, choose from min_len to end_pos
            if overhang not in overhangs_1st or not check_overhang_unique(overhang_list, overhang):
                frag, rest_seq, overhang, overhang_pair = find_overhang_from_list_min_re(end_pos, frag_min_len, overhang_list, seq, overhang_len)

            if overhang not in overhangs_1st or not check_overhang_unique(overhang_list, overhang):
                frag, rest_seq, overhang, overhang_pair = find_overhang_from_2nd_re(end_pos, i, max_len_3, max_inner_len, frag_min_len, overhang_list, seq, overhang_len, frag_num)
            
        ##if overhang in overhangs_1st:
        #     print("choose from best group: ", overhang)
        # else:
        #     print("overhang:", overhang )
        # if not check_overhang_unique(overhang_list, overhang):
        #     print("no unique")
    if split_successful_tag:
        return frag, rest_seq, overhang, overhang_pair
    else:
        return split_successful_tag
    
def add_adapter_list(seq_name, frag_seq, adapter_list, adapter_idx, frag_num, adapter_F, adapter_R, seq_barcode):
    # It can add adapter to frags
    alphabet = ["a", "b", "c", "d", "e", "f", "g"]
    frag_ad_list = {}
    for i in range(frag_num):
        name = alphabet[i]+ "_"+ seq_name
        if frag_seq[i] == "":
            frag_ad_list[name] = "not_split"
        else:
            if i == 0:
                frag_ad = adapter_F + adapter_list[i][adapter_idx] + NcoI + seq_barcode[0] + NdeI + frag_seq[i] + adapter_list[i + 1][adapter_idx] + adapter_R
                
                if len(frag_ad) < 211 :
                    frag_ad = adapter_F + adapter_list[i][adapter_idx] + NcoI + seq_barcode[0] + NdeI + frag_seq[i] + 'AAAAAAA' + 'agggctcgtgtcgccaccaatggggtattcacgtagcggctgggctcgtgtcgccaccaatggggtattcacgtagcggctgg' + adapter_list[i + 1][adapter_idx] + adapter_R
                elif len(frag_ad) < 251 :
                    frag_ad = adapter_F + adapter_list[i][adapter_idx] + NcoI + seq_barcode[0] + NdeI + frag_seq[i] + 'AAAAAAA' + 'agggctcgtgtcgccaccaatggggtattcacgtagcggctgg' + adapter_list[i + 1][adapter_idx] + adapter_R
                
                frag_ad_list[name] = frag_ad
                
            elif i == frag_num - 1:
                frag_ad = adapter_F + adapter_list[i][adapter_idx] + frag_seq[i] + XhoI + seq_barcode[1] + BamHI + adapter_list[i + 1][adapter_idx] + adapter_R
                
                if len(frag_ad) <211:
                    frag_ad = adapter_F + adapter_list[i][adapter_idx] + 'agggctcgtgtcgccaccaatggggtattcacgtagcggctgggctcgtgtcgccaccaatggggtattcacgtagcggctgg' + 'CCCCCCC' + frag_seq[i] + XhoI + seq_barcode[1] + BamHI + adapter_list[i + 1][adapter_idx] + adapter_R
                elif len(frag_ad) <251:
                    frag_ad = adapter_F + adapter_list[i][adapter_idx] + 'agggctcgtgtcgccaccaatggggtattcacgtagcggctgg' + 'CCCCCCC' + frag_seq[i] + XhoI + seq_barcode[1] + BamHI + adapter_list[i + 1][adapter_idx] + adapter_R
                
                frag_ad_list[name] = frag_ad
            
            else:
                frag_ad = adapter_F + adapter_list[i][adapter_idx] + frag_seq[i] + adapter_list[i + 1][adapter_idx] + adapter_R
                
                if len(frag_ad) <211:
                    frag_ad = adapter_F + adapter_list[i][adapter_idx] + frag_seq[i] + "enzyme" + 'agggctcgtgtcgccaccaatggggtattcacgtagcggctgggctcgtgtcgccaccaatggggtattcacgtagcggctgg' + adapter_list[i + 1][adapter_idx] + adapter_R
                elif len(frag_ad) < 251 :
                    frag_ad = adapter_F + adapter_list[i][adapter_idx] + frag_seq[i] + "enzyme" + 'agggctcgtgtcgccaccaatggggtattcacgtagcggctgg' + adapter_list[i + 1][adapter_idx] + adapter_R
                
                frag_ad_list[name] = frag_ad
 
    return frag_ad_list

def split_sequences(seq_list, frag_num, adapter_list, adapter_idx, adapter_F, adapter_R, seq_barcode_list, max_oligo_size=300, overhang_len=4):
    assert len(seq_list) < 120
    
    primer_len_5=len(adapter_list[0][adapter_idx])# 5端primer长度
    primer_len_3=len(adapter_list[-1][adapter_idx])# 3端primer长度
    primer_in_len=len(adapter_list[1][adapter_idx])# 中间primer长度 基于a片段的3端
    adapter_F_len = len(adapter_F) #5'端adapter长度
    adapter_R_len = len(adapter_R) #3'端adapter长度
    seq_barc_len = len(seq_barcode_list[0][0]) #sequence barcode length

    max_len_5 = max_oligo_size - adapter_F_len - 6 - primer_len_5 - seq_barc_len  - 6 - primer_in_len  - adapter_R_len
    #5端frag最长长度
    max_inner_len = max_oligo_size - adapter_F_len - primer_in_len * 2 - adapter_R_len
    #中间frag最长长度
    max_len_3 = max_oligo_size - adapter_F_len - primer_in_len - 6 - seq_barc_len - primer_len_3 - 6 - adapter_R_len
    #3端frag最长长度
    
    cut_successful_flag = False
    overhang_list = []
    while not cut_successful_flag:
        for frag_order in range(1, frag_num):
            frag_order = []
            overhang_list.append(frag_order)

        split_result = {}
        
        seq_num = 0
        for name, seq in tqdm(seq_list.items()):
            not_unique_overhang = False
            frag_seq = []
            if seq_num % 2 == 0:
                for i in range(1, frag_num):
                    frag_num_left = frag_num - i
                    frag_min_len=count_frag_lenth(len(seq), i, frag_num_left, max_inner_len=max_inner_len, max_fist_len=max_len_5, max_last_len=max_len_3)
                    
                    #End_pos is the average length of each fragment.
                    end_pos = round((len(seq) + frag_num_left * 4 )/(frag_num_left+1))
                    
                    result = find_unique_overhang(end_pos, i, max_len_5, max_len_3, max_inner_len, frag_min_len, overhang_list[i-1], seq, seq_num, frag_num)
                    #find unique overhang
                    if isinstance(result, bool):
                        cut_successful_flag = True
                    else:
                        frag, rest_seq, overhang, overhang_pair = result
                    #Put frag into a list 
                    while not check_overhang_unique(overhang_list[i-1], overhang):
                        not_unique_overhang = True
                        frag_seq.append("no")
                        break
                    else:
                        frag_seq.append(frag)
                    if i == frag_num - 1:
                        frag_seq.append(rest_seq)
                    if not_unique_overhang:
                        break

                    overhang_list[i-1].append(overhang)
                    overhang_list[i-1].append(overhang_pair)
                    seq = rest_seq 
            else:
                for i in range(frag_num - 1, 0, -1):
                    frag_num_left = i
                    frag_order = frag_num - i
                    frag_min_len=count_frag_lenth(len(seq), frag_order, frag_num_left, max_inner_len=max_inner_len, max_fist_len=max_len_5, max_last_len=max_len_3)
                    
                    #End_pos is the average length of each fragment.
                    aver_lenth = round((len(seq) + frag_num_left * 4 )/(frag_num_left+1))
                    end_pos = len(seq) - aver_lenth

                    #find unique overhang
                    frag, rest_seq, overhang, overhang_pair = find_unique_overhang(end_pos, i, max_len_5, max_len_3, max_inner_len, frag_min_len, overhang_list[i-1], seq, seq_num, frag_num)
                    #Put frag into a list 
                    while not check_overhang_unique(overhang_list[i-1], overhang):
                        print("no unique")
                        not_unique_overhang = True
                        frag_seq.insert(0, "no")
                        break
                    else:
                        frag_seq.insert(0, frag)
                    if i == 1:
                        frag_seq.insert(0, rest_seq)
                    if not_unique_overhang:
                        break

                    overhang_list[i-1].append(overhang)
                    overhang_list[i-1].append(overhang_pair)
                    seq = rest_seq 
            
            no_split = True if ("no" in frag_seq) else print("Oligo %s has been successfully spilted." %(name))
            if not no_split:
                seq_barcode = seq_barcode_list[seq_num]
                frag_list = add_adapter_list(name, frag_seq, adapter_list, adapter_idx, frag_num, adapter_F, adapter_R, seq_barcode) 
                split_result[name]=frag_list 
            else:
                print("Oligo %s has not been successfully spilted." %(name))
            
            seq_num += 1
        if not_unique_overhang:
            continue
        else:
            print("Work is done!\n")
            cut_successful_flag = True
    while len(split_result) < len(seq_list):
        return False
    else:
        return split_result, overhang_list

def main(args):
    ####################################
    #input oligo length
    max_oligo_size = args.max_oligo_size
    min_len = args.min_oligo_size
    frag_num = args.frag_num
    ad_fname = args.adapter_fname
    seq_barc_fname = args.sequence_barcode_fname
    adapter_idx = args.adapter_number - 1
    adapter_F = args.adapter_f
    adapter_R = args.adapter_r

    #input the design names and sequences
    with open(args.input_list,'r', encoding='utf-8-sig') as ip:
        lines = ip.readlines()

    seq_list = {}
    for line in lines:
        its = line.strip().split()
        if len(its) != 2:
            #print("Skipping:")
            #print(its)
            continue
        else:
            seq_list[its[0]] = its[1]

    #print(f'file loaded! {len(seq_list)}')

    '''
    BEGGINING OF MAIN:
    '''
    # Load and process sequence files:
    # for this verison, specify which line of the adapters should be used
    # numbering starts from 1, not like python 0
    
    # Get adapters ready
    adapter_list = read_adapters(ad_fname, frag_num)
    
    # Get sequence barcodes ready
    seq_barcode_list = read_sequence_barcode(seq_barc_fname, len(seq_list))
    
    #split the sequences into oligos
    result = split_sequences(seq_list, frag_num, adapter_list, adapter_idx, adapter_F, adapter_R, seq_barcode_list)
    while isinstance(result, bool):
        seq_list = shuffle_dict(seq_list)
        result = split_sequences(seq_list, frag_num, adapter_list, adapter_idx, adapter_F, adapter_R, seq_barcode_list)
    else:
        result_list, overhang_list = result
        
    #write the output file
    with open('%d_oligos_subpools_%d.tab'%(frag_num, args.adapter_number), 'w') as outputfile:
        for name, frag_seq in result_list.items():
            for key_1, value_1 in frag_seq.items():
                if len(value_1) > 300:
                    print('Oligo %s is longer than 300bp!!!' % key_1)
                if len(value_1) < 250:
                    print('Oligo %s is shorter than 250bp!!!' % key_1)
                outputfile.write(key_1 + ', ' + value_1 + '\n')


if __name__ == '__main__':
    ########## Option system: ###########
    argparser = argparse.ArgumentParser(description='Split genes in orthogonal pieces that can be used in multiplex assembly')
    argparser.add_argument('-input_list', type=str, help='Name of file containing: mygenename DNAsequence ')
    argparser.add_argument('-adapter_fname', type=str, default='./pool_adapters_short_list.txt',help='Name of file containing adapter sequences. Format: First line: column names, followed by lines: adapter_name fiveprime_5 fiveprime_3 threeprime_5 threeprime_3')
    argparser.add_argument('-adapter_f', type=str, default='FFFFFFFFFFFFFFFFFFFF',help='Forward adapter sequences. Default: ')
    argparser.add_argument('-adapter_r', type=str, default='RRRRRRRRRRRRRRRRRRRR',help='Reverse adapter sequences. Default: ')
    argparser.add_argument('-sequence_barcode_fname', type=str, default='./sequence_barcode_list.txt',help='Name of file containing sequence barcode sequences. Format: First line: column names, followed by lines: fiveprime_5 threeprime_3')
    argparser.add_argument('-adapter_number', type=int, help='What adapter to use? starting at 1')
    argparser.add_argument('-max_oligo_size', type=int, default=300, help='Absolute max length of orderable oligo')
    argparser.add_argument('-min_oligo_size', type=int, default=120, help='Absolute min length of orderable oligo')
    argparser.add_argument('-enzyme_num', type=int, default=1, choices=[1, 2], help='The number of how many enzymes.')
    argparser.add_argument('-frag_num', type=int, default=3, help='The number of how many fragments you want to split.')
    
    args = argparser.parse_args()
    main(args)

