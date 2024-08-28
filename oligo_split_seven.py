#!/usr/bin/env python
# Mar 28, 2024
# Zhien Wu modified
import numpy as np
from tqdm import tqdm
import sys
import random
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

#overhand_list
#The best overhangs are included in overhangs_1st. The second good overhangs are included in bad_overhangs.
overhangs_1st = ['TCAA', 'ATAA', 'TAGA', 'GTTA', 'AGGG', 'CCTA', 'AAGA', 'TCCA', 'AGAA', 'AAAG', 'ACAT', 'GGGA', 'GCAA', 'TAAA', 'TGAA', 'GACA', 'ACGA', 'GGTA', 'ATCC', 'ATTG', 'CTAC', 'AAGC', 'CATC', 'ACTC', 'CACA', 'CTAA', 'GAAA', 'AGAC', 'AGCA', 'CGTA', 'ACCG', 'CAGA', 'AACT', 'AATA', 'GCAC', 'CCAG', 'CAAG', 'AAAT', 'ATCA', 'CAGG', 'CATA', 'GGAA', 'AGGA', 'ACGC', 'ATAC', 'CTCA', 'GCCA', 'CCGA', 'ACAG', 'AATC', 'CAGC', 'AAAA', 'AGTG', 'CGCA', 'AACG', 'GAGA', 'ACTA', 'TACA', 'ATGA', 'CGAC', 'CGAA', 'AGCC']
bad_overhangs = ['GGGG', 'CCCC', 'GGGC', 'GCCC', 'GGCG', 'CGCC', 'GCGG', 'CCGC', 'CGGG', 'CCCG', 'GGCC', 'CCGG', 'GCCG', 'CGGC', 'CGCG', 'GCGC', 'ATTA', 'TAAT']#too much?

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

#function:
def shuffle_dict(input_dict):
    keys = list(input_dict.keys())
    random.shuffle(keys)
    shuffled_dict = {key: input_dict[key] for key in keys}
    return shuffled_dict

def read_subpool_barcode(infile, frag_num):
    """
    Read the subpool barcode file, obetaining the correct sequences from the expected format.
    """
    handle = open(infile, 'r')
    # get all lines, skip final newline and drop first line that has titles:
    lines = [ line[:-1] for line in handle.readlines() ][1:]
    spool_barcode_list=[]
    for i in range(frag_num * 2):
        l = [ line.split()[i+1] for line in lines ]
        spool_barcode_list.append(l)
    return spool_barcode_list

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

def count_frag_lenth(seq_len, frag_th, frag_num_left, max_inner_len=220, max_last_len=210, max_fist_len=210, min_len=115):
    #This function is used to calculate the minimum length of the fragment.
    
    frag_min_len = seq_len - (max_inner_len - 4) * ( frag_num_left - 1 ) - ( max_last_len - 4)
    
    if frag_min_len < min_len:
        frag_min_len = min_len
    
    if frag_th == 1:
        assert frag_min_len <= max_fist_len
    else:
        assert frag_min_len <= max_inner_len

    if frag_num_left == 1:  #haven't used
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

def find_overhang_from_list_max(end_pos, i, max_len_5, max_inner_len, overhang_list, seq, overhang_len):
    reasonable_tag = True
    
    max_len = max_len_5 if i == 1 else max_inner_len
    if max_len < end_pos:
        reasonable_tag = False
        return reasonable_tag
    else:
        for j in range(end_pos, max_len + 1):
            frag, rest_seq, overhang, overhang_pair = split_with_pos(seq, j, overhang_len)
            if overhang in overhangs_1st and check_overhang_unique(overhang_list, overhang):
                break
        
        #If there is no unique overhang from end_pos to max_len, return the last overhang
        return frag, rest_seq, overhang, overhang_pair
   
def find_overhang_from_list_min(end_pos, frag_min_len, overhang_list, seq, overhang_len):
    reasonable_tag = True
    
    if frag_min_len > end_pos:
        reasonable_tag = False
        return reasonable_tag
    else:
        for j in range(end_pos, frag_min_len - 1, -1):
            frag, rest_seq, overhang, overhang_pair = split_with_pos(seq, j, overhang_len)
            if overhang in overhangs_1st and check_overhang_unique(overhang_list, overhang):
                break
        
        #If there is no unique overhang from end_pos to max_len, return the last overhang
        return frag, rest_seq, overhang, overhang_pair

def find_overhang_from_2nd(end_pos, i, max_len_5, max_inner_len, frag_min_len, overhang_list, seq, overhang_len):
    reasonable_tag = True
    
    max_len = max_len_5 if i == 1 else max_inner_len
    if max_len < end_pos:
        reasonable_tag = False
        return reasonable_tag
    else:
        for j in range(end_pos, max_len + 1):
            frag, rest_seq, overhang, overhang_pair = split_with_pos(seq, j, overhang_len)
            if overhang not in bad_overhangs and check_overhang_unique(overhang_list, overhang):
                return frag, rest_seq, overhang, overhang_pair
                break
            #If there is no unique overhang from end_pos to max_len, overhang is the last overhang
            
        while overhang in bad_overhangs or not check_overhang_unique(overhang_list, overhang):
            if end_pos < frag_min_len:
                reasonable_tag = False
                return reasonable_tag
                break
            for j in range(frag_min_len, end_pos):
                frag, rest_seq, overhang, overhang_pair = split_with_pos(seq, j, overhang_len)
                if overhang not in bad_overhangs and check_overhang_unique(overhang_list, overhang):
                    return frag, rest_seq, overhang, overhang_pair
                    break
                
                #If there is no unique overhang, return the resonable tag(True)
                return reasonable_tag

def find_unique_overhang(end_pos, 
                         i, 
                         max_len_5, 
                         max_len_3, 
                         max_inner_len, 
                         frag_min_len, 
                         overhang_list, 
                         seq, 
                         seq_num, 
                         frag_num, 
                         random_choice = None, 
                         overhang_len = 4):
    
    split_successful_tag = False
    
    while not split_successful_tag:
        #Randomly choose to start from max_len or min_len. If random_choice is 0, start from max_len. If random_choice is 1, start from min_len. While random_choice has been chosen, it will be changed into the other one.
        if random_choice == None:
            random_choice = np.random.randint(0,2)
        elif random_choice == 0:
            random_choice = 1
        else:
            random_choice = 0
            
        if random_choice == 0: 
            #Choose overhang from end_pos to max_len
            result = find_overhang_from_list_max(end_pos, i, max_len_5, max_inner_len, overhang_list, seq, overhang_len)
            if isinstance(result, bool):
                break
            else:
                frag, rest_seq, overhang, overhang_pair = result
                
            if overhang in overhangs_1st and check_overhang_unique(overhang_list, overhang):
                split_successful_tag = True
            else:
                #If there is no unique overhang from end_pos to max_len, choose from min_len to end_pos
                result = find_overhang_from_list_min(end_pos, frag_min_len, overhang_list, seq, overhang_len)
                if isinstance(result, bool):
                    break
                else:
                    frag, rest_seq, overhang, overhang_pair = result
                    
                #If there is no unique overhang from 1st list, choose from min_len to max_len      
                if overhang in overhangs_1st and check_overhang_unique(overhang_list, overhang):
                    split_successful_tag = True
                else:
                    result = find_overhang_from_2nd(end_pos, i, max_len_5, max_inner_len, frag_min_len, overhang_list, seq, overhang_len)
                    if isinstance(result, bool):
                        break
                    else:
                        frag, rest_seq, overhang, overhang_pair = result
                        split_successful_tag = True
                
        else:
            #Choose overhang from min_len to end_pos
            result = find_overhang_from_list_min(end_pos, frag_min_len, overhang_list, seq, overhang_len)
            if isinstance(result, bool):
                break
            else:
                frag, rest_seq, overhang, overhang_pair = result
                
            if overhang in overhangs_1st and check_overhang_unique(overhang_list, overhang):
                split_successful_tag = True
            else:
            #If there is no unique overhang from min_len to end_pos, choose from end_pos to max_len
                result = find_overhang_from_list_max(end_pos, i, max_len_5, max_inner_len, overhang_list, seq, overhang_len)
                if isinstance(result, bool):
                    break
                else:
                    frag, rest_seq, overhang, overhang_pair = result
                    
                if overhang in overhangs_1st and check_overhang_unique(overhang_list, overhang):
                    split_successful_tag = True
                else:
                    #If there is no unique overhang from 1st list, choose from min_len to max_len
                    result = find_overhang_from_2nd(end_pos, i, max_len_5, max_inner_len, frag_min_len, overhang_list, seq, overhang_len)
                    if isinstance(result, bool):
                        break
                    else:
                        frag, rest_seq, overhang, overhang_pair = result
                        if check_overhang_unique(overhang_list, overhang):
                            split_successful_tag = True
                        else:
                            split_successful_tag = False
            
    if split_successful_tag:
        return frag, rest_seq, overhang, overhang_pair, random_choice
    else:
        return split_successful_tag
    
def add_spool_barcode_list(seq_name, 
                           frag_seq, 
                           spool_barcode_list, 
                           subp_barc_idx, 
                           frag_num, 
                           adapter_F, 
                           adapter_R, 
                           seq_barcode):
    
    # It can add adapter to frags
    alphabet = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z"]
    frag_ad_list = {}
    
    for i in range(frag_num):
        name = alphabet[i]+ "_"+ seq_name
        if frag_seq[i] == "":
            frag_ad_list[name] = "not_split"
        else:
            if i == 0:
                frag_ad = adapter_F + spool_barcode_list[2 * i][subp_barc_idx] + NcoI + seq_barcode[0] + NdeI + frag_seq[i] + bsai + spool_barcode_list[2*i + 1][subp_barc_idx] + adapter_R
                
                if len(frag_ad) < 211 :
                    frag_ad = adapter_F + spool_barcode_list[2 * i][subp_barc_idx] + NcoI + seq_barcode[0] + NdeI + frag_seq[i] + bsai + 'atgagccatattcaacgggaaacgtcttgctgcgattaaattccaacatggatgctgatttatatgggtatataat' + bsai + spool_barcode_list[2*i + 1][subp_barc_idx] + adapter_R
                elif len(frag_ad) < 251 :
                    frag_ad = adapter_F + spool_barcode_list[2 * i][subp_barc_idx] + NcoI + seq_barcode[0] + NdeI + frag_seq[i] + bsai + 'atgagccatattcaacgggaaacgtcttgctgtaat' + bsai + spool_barcode_list[2*i + 1][subp_barc_idx] + adapter_R
                
                frag_ad_list[name] = frag_ad
                
            elif i == frag_num - 1:
                #the last fragment
                if i % 2 == 0:
                    enzyme = BsmBI
                else:
                    enzyme = BsaI
                    
                frag_ad = adapter_F + spool_barcode_list[2 * i][subp_barc_idx] + enzyme + frag_seq[i] + XhoI + seq_barcode[1] + BamHI + spool_barcode_list[2*i + 1][subp_barc_idx] + adapter_R
                
                if len(frag_ad) < 211:
                    frag_ad = adapter_F + spool_barcode_list[2 * i][subp_barc_idx] + enzyme + 'atgagccatattcaacgggaaacgtcttgctgcgattaaattccaacatggatgctgatttatatgggtatataat' + enzyme + frag_seq[i] + XhoI + seq_barcode[1] + BamHI + spool_barcode_list[2*i + 1][subp_barc_idx] + adapter_R
                elif len(frag_ad) <251:
                    frag_ad = adapter_F + spool_barcode_list[2 * i][subp_barc_idx] + enzyme + 'atgagccatattcaacgggaaacgtcttgctgtaat' + enzyme + frag_seq[i] + XhoI + seq_barcode[1] + BamHI + spool_barcode_list[2*i + 1][subp_barc_idx] + adapter_R
                
                frag_ad_list[name] = frag_ad
            
            else:
                if i % 2 == 0:
                    enzyme_5 = BsmBI
                    enzyme_3 = bsai
                else:
                    enzyme_5 = BsaI
                    enzyme_3 = bsmbi
                
                frag_ad = adapter_F + spool_barcode_list[2 * i][subp_barc_idx] + enzyme_5 + frag_seq[i] + enzyme_3 + spool_barcode_list[2*i + 1][subp_barc_idx] + adapter_R
                
                if len(frag_ad) < 211:
                    frag_ad = adapter_F + spool_barcode_list[2 * i][subp_barc_idx] + enzyme_5 + frag_seq[i] + enzyme_3 + 'atgagccatattcaacgggaaacgtcttgctgcgattaaattccaacatggatgctgatttatatgggtatataat' + enzyme_3 + spool_barcode_list[2*i + 1][subp_barc_idx] + adapter_R
                elif len(frag_ad) < 251 :
                    frag_ad = adapter_F + spool_barcode_list[2 * i][subp_barc_idx] + enzyme_5 + frag_seq[i] + enzyme_3 + 'atgagccatattcaacgggaaacgtcttgctgtaat' + enzyme_3 + spool_barcode_list[2*i + 1][subp_barc_idx] + adapter_R
                
                frag_ad_list[name] = frag_ad
 
    return frag_ad_list

def split_sequences(seq_list, 
                    frag_num, 
                    spool_barcode_list, 
                    subp_barc_idx, 
                    adapter_F, 
                    adapter_R, 
                    seq_barcode_list, 
                    max_oligo_size=301,
                    min_oligo_len=251):
    
    #Check if the size of subpool is larger than 100
    assert len(seq_list) < 100
    
    #subpool barcode length, all the same
    subp_barc_len_5 = len( spool_barcode_list[0][subp_barc_idx] ) #  the length of 5' subpool barcode
    subp_barc_len_3 = len( spool_barcode_list[-1][subp_barc_idx] ) #  the length of 3' subpool barcode
    subp_barc_in_len = len( spool_barcode_list[1][subp_barc_idx] ) #  the length of inner subpool barcode
    #adapter length
    adapter_F_len = len( adapter_F ) #  the length of 5' adapter
    adapter_R_len = len( adapter_R ) #  the length of 3' adapter
    #sequence barcode length
    seq_barc_len = len(seq_barcode_list[0][0]) #  the length of sequence barcode

    #the length of the longest A fragment
    max_len_5 = max_oligo_size - adapter_F_len - subp_barc_len_5 - 6 - seq_barc_len  - 6 - (subp_barc_in_len + 7) - adapter_R_len #  add the basi or not?
    #the length of the longest inner fragment
    max_inner_len = max_oligo_size - adapter_F_len - (subp_barc_in_len + 7) * 2 - adapter_R_len
    #the length of the longest last fragment
    max_len_3 = max_oligo_size - adapter_F_len - (subp_barc_in_len + 7) - 6 - seq_barc_len - 6 - subp_barc_len_3 - adapter_R_len
    
    #generate overhang list
    overhang_list = []
    for frag_order in range(1, frag_num):
        frag_order = []
        overhang_list.append(frag_order)

    split_result = {}
    
    #split the sequences
    seq_num = 0
    for name, seq in tqdm(seq_list.items()):
        
        #keep the original sequence
        init_seq = seq
        cut_successful_flag = False
        not_unique_overhang = False
        
        try_times = 0
        while not cut_successful_flag:
            frag_seq = []
            seq = init_seq
            for i in range(1, frag_num):
                frag_num_left = frag_num - i
                frag_min_len = count_frag_lenth(len(seq), 
                                                i, 
                                                frag_num_left, 
                                                max_inner_len=max_inner_len, 
                                                max_fist_len=max_len_5, 
                                                max_last_len=max_len_3)
                
                #End_pos is the average length of each fragment.
                end_pos = round((len(seq) + frag_num_left * 4 )/(frag_num_left+1))

                #find unique overhang
                result = find_unique_overhang(end_pos, 
                                              i, 
                                              max_len_5, 
                                              max_len_3, 
                                              max_inner_len, 
                                              frag_min_len, 
                                              overhang_list[i-1], 
                                              seq, 
                                              seq_num, 
                                              frag_num)
                
                if isinstance(result, bool):
                    #can't find unique overhang, need to re-split
                    try_times += 1
                    break
                else:
                    frag, rest_seq, overhang, overhang_pair, random_choice = result
                    #Put frag into a list 
                    frag_seq.append(frag)
                    if i == frag_num - 1:
                        frag_seq.append(rest_seq)
                    
                    #Add overhang to the list
                    overhang_list[i-1].append(overhang)
                    overhang_list[i-1].append(overhang_pair)
                    seq = rest_seq 

            #If the sequence can't be split after 5 times, skip this sequence
            if isinstance(result, bool) and try_times < 5:
                continue
            if try_times == 5:
                not_unique_overhang = True
            else:
                cut_successful_flag = True
                
        if not_unique_overhang:
            print("Oligo %s" %(name))
        else:
            seq_barcode = seq_barcode_list[seq_num]
            frag_list = add_spool_barcode_list(name, 
                                               frag_seq, 
                                               spool_barcode_list, 
                                               subp_barc_idx, 
                                               frag_num, 
                                               adapter_F, 
                                               adapter_R, 
                                               seq_barcode) 
            split_result[name]=frag_list
        
        seq_num += 1
    
    #Check if all sequences have length between 250 and 300
    not_suitable_length_flag = False
    for name, frag_seq in split_result.items():
        for key_1, value_1 in frag_seq.items():
            if len(value_1) < min_oligo_len or len(value_1) > max_oligo_size:
                not_suitable_length_flag = True
    #Check if all sequences have been split
    if len(split_result) < len(seq_list) or not_suitable_length_flag:
        print("Try again!")
        return False
    else:
        return split_result, overhang_list

def main(args):
    ####################################
    #input oligo length
    min_oligo_len = args.min_oligo_length
    frag_num = args.frag_num
    spool_barc_fname = args.subpool_barcode_fname
    seq_barc_fname = args.sequence_barcode_fname
    subp_barc_idx = args.subp_barc_index - 1
    adapter_F = args.adapter_f
    adapter_R = args.adapter_r

    #input the design names and sequences
    with open(args.input_list,'r', encoding='utf-8-sig') as ip:
        lines = ip.readlines()

    seq_list = {}
    for line in lines:
        its = line.strip().split()
        if len(its) != 2:
            continue
        else:
            seq_list[its[0]] = its[1]

    '''
    BEGGINING OF MAIN:
    '''
    # Load and process sequence files:
    # for this verison, specify which line of the adapters should be used
    # numbering starts from 1, not like python 0
    
    # Get adapters ready
    spool_barcode_list = read_subpool_barcode(spool_barc_fname, frag_num)
    
    # Get sequence barcodes ready
    seq_barcode_list = read_sequence_barcode(seq_barc_fname, len(seq_list))
    
    #split the sequences into oligos
    result = split_sequences(seq_list, 
                             frag_num, 
                             spool_barcode_list, 
                             subp_barc_idx, 
                             adapter_F, adapter_R, 
                             seq_barcode_list, 
                             max_oligo_size = args.max_oligo_length,
                             min_oligo_len = args.min_oligo_length)
    
    while isinstance(result, bool):
        seq_list = shuffle_dict(seq_list)
        result = split_sequences(seq_list, 
                                 frag_num, 
                                 spool_barcode_list, 
                                 subp_barc_idx, 
                                 adapter_F, adapter_R, 
                                 seq_barcode_list, 
                                 max_oligo_size = args.max_oligo_length,
                                 min_oligo_len = args.min_oligo_length)
    else:
        result_list, overhang_list = result
    
    #write the output file
    with open('%d_oligos_subpools_%d.tab'%(frag_num, args.subp_barc_index), 'w') as outputfile:
        for name, frag_seq in result_list.items():
            for key_1, value_1 in frag_seq.items():
                if len(value_1) > args.max_oligo_length:
                   print('Oligo %s is longer than %dbp!!!' % (key_1, args.max_oligo_length))
                if len(value_1) < args.min_oligo_length:
                    print('Oligo %s is shorter than %dbp!!!' % (key_1, args.min_oligo_length))
                outputfile.write(key_1 + ', ' + value_1 + '\n')


if __name__ == '__main__':
    ########## Option system: ###########
    argparser = argparse.ArgumentParser(description='Split genes in orthogonal pieces that can be used in multiplex assembly')
    argparser.add_argument('--input_list', type=str, help='Name of file containing: mygenename DNAsequence ')
    argparser.add_argument('--subpool_barcode_fname', type=str, default='./pool_subpools_barcode.txt',help='Name of file containing adapter sequences. Format: First line: column names, followed by lines: adapter_name fiveprime_5 fiveprime_3 threeprime_5 threeprime_3')
    argparser.add_argument('--adapter_f', type=str, default='FFFFFFFFFFFFFFFFFFFF',help='Forward adapter sequences. Default: ')
    argparser.add_argument('--adapter_r', type=str, default='RRRRRRRRRRRRRRRRRRRR',help='Reverse adapter sequences. Default: ')
    argparser.add_argument('--sequence_barcode_fname', type=str, default='./sequence_barcode_list.txt',help='Name of file containing sequence barcode sequences. Format: First line: column names, followed by lines: fiveprime_5 threeprime_3')
    argparser.add_argument('--subp_barc_index', type=int, help='What subpool barcode to use? starting at 1')
    argparser.add_argument('--max_oligo_length', type=int, default=300, help='Absolute max length of orderable oligo')
    argparser.add_argument('--min_oligo_length', type=int, default=251, help='Absolute min length of orderable oligo')
    argparser.add_argument('--enzyme_num', type=int, default=1, choices=[1, 2], help='The number of how many enzymes.')
    argparser.add_argument('--frag_num', type=int, default=3, help='The number of how many fragments you want to split.')
    
    args = argparser.parse_args()
    main(args)

