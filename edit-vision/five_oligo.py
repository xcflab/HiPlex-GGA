#!/usr/bin/env python
import numpy as np
from tqdm import tqdm
import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

#function:
def read_adapters(infile, frag_num):
    """
    Reads in the adapter file, obetaining the correct sequences from the expected format.
    """
    handle = open(infile, 'r')
    # get all lines, skip final newline and drop first line that has titles:
    lines = [ line[:-1] for line in handle.readlines() ][1:]
    # print("Your input adapters:")
    # for line in lines:
    #     print(line)
    adapter_list=[]
    for i in range(frag_num * 2):
        l = [ line.split()[i+1] for line in lines ]
        adapter_list.append(l)
    #oligoa_5 = [ line.split()[1] for line in lines ]
    #oligoa_3 = [ line.split()[2] for line in lines ]
    #oligob_5 = [ line.split()[3] for line in lines ]
    #oligob_3 = [ line.split()[4] for line in lines ]
    #oligoc_5 = [ line.split()[5] for line in lines ]
    #oligoc_3 = [ line.split()[-1] for line in lines ]
    #return oligoa_5, oligoa_3, oligob_5, oligob_3, oligoc_5, oligoc_3
    return adapter_list

def count_frag_lenth(seq_len, frag_th, max_inner_len, frag_num_left, max_last_len, max_fist_len, min_len=150):
    print("total length:", seq_len)
    frag_min_len = seq_len - (max_inner_len - 4) * ( frag_num_left - 1 ) - ( max_last_len - 4)
    print("frag min length:", frag_min_len)
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
    # 根据粘性末端的位置切分序列
    end_len_1 = overhang_len
    seq_1 = sequence[0:pos]
    seq_2 = sequence[pos-end_len_1:]
    overhang = Seq(sequence[pos-end_len_1:pos])
    overhang_pair = overhang.reverse_complement()
    return seq_1, seq_2, overhang, overhang_pair

def check_ruslt(frag_1,frag_2,frag_3):
    # 用来检查结果的可视化
    f1_end = frag_1[-10:]
    f2_start = frag_2[:10]
    f2_end = frag_2[-10:]
    f3_start = frag_3[:10]

    bar = "-"*6
    print("\ncohesive end 1")
    print(f"{f1_end}{bar}")
    print(f"{bar}{f2_start}")

    print("\ncohesive end 2")
    print(f"{f2_end}{bar}")
    print(f"{bar}{f3_start}")
    
def add_adapter(frag, adapter_5, adapter_3, enzyme):
    # It can add adapter to frags
    frag_ada = adapter_5 + frag + adapter_3
    if len(frag_ada) < 201 :
        frag_ada = adapter_5 + frag + enzyme + 'cactgcgtggttcgcgtcctaaaccagtggccgggatagacacttgatccaaattgtgactcaccaaaggtacatcgtcacattgccaactgg' + adapter_3
    elif len(frag_ada) < 251 :
        frag_ada = adapter_5 + frag + enzyme + 'agggctcgtgtcgccaccaatggggtattcacgtagcggctgg' + adapter_3
 
    return frag_ada

def add_adapter_last(frag, adapter_5, adapter_3, enzyme):
    # It can add adapter to frags
    frag_ada = adapter_5 + frag + adapter_3
    if len(frag_ada) < 201 :
        frag_ada = adapter_5 + 'ctggctctaagtaggaaccagacgagatacgttattgactttaactccccagcagtttgtgagccttgcgcggtgtccgaagactctcgacga' + enzyme + frag + adapter_3
    elif len(frag_ada) < 251 :
        frag_ada = adapter_5 + 'ctggttttaatcagacgccacctctctcgcgagacacaatagg' + enzyme + frag + adapter_3
    return frag_ada

def split_sequences_more_fragments(seq_list, frag_num, adapter_list, adapter_idx, max_oligo_size=300, primer_len_5=50, primer_in_len=40, primer_len_3=50, overhang_len=4):
    '''
    使用两种酶的情况下, 根据oligo长度和分段切分末端
    '''
    #input the site of enzyme
    #3'
    bsai = 'AGAGACC'
    bsmi = 'AGAGACG'
    #5'
    BsaI = 'GGTCTCA'
    BsmI = 'CGTCTCA'
    
    oligoa_3 = adapter_list[1]
    oligob_5 = adapter_list[2]
    oligob_3 = adapter_list[3]
    oligoc_5 = adapter_list[4]
    ecoli_primer5 = adapter_list[0][adapter_idx]# 5' == adapter_list[0]
    ecoli_primer3 = adapter_list[-1][adapter_idx] # 3' == adapter_list[-1]
    primer_len_5=len(ecoli_primer5)# 5端primer长度
    primer_len_3=len(ecoli_primer3)# 3端primer长度
    primer_in_len=len(oligoa_3[adapter_idx])# 中间primer长度 基于a片段的3端
    batch_big_id = adapter_idx

    max_len_5 = ( max_oligo_size - primer_len_5 - primer_in_len )
    #5端frag最长长度
    max_inner_len =( max_oligo_size - primer_in_len * 2 )
    #中间frag最长长度
    max_len_3 = ( max_oligo_size - primer_len_3 - primer_in_len )
    #3端frag最长长度
    #print(max_len_5, max_inner_len, max_len_3)
        
    assert len(seq_list) < 101
    cut_successful_flag = False
    overhang_list = []
    while not cut_successful_flag:
        for frag_order in range(1, frag_num):
            frag_order = [Seq("CTGG")]
            overhang_list.append(frag_order)
        
        split_result = {}

        for name, seq in tqdm(seq_list.items()):
            print(name)
            for i in range(1, frag_num):
                print(i)
                frag_num_left = frag_num - i

                frag_min_len=count_frag_lenth(len(seq), i, max_inner_len, frag_num_left, max_last_len=max_len_3, max_fist_len=max_len_5)
                
                # 初始化为n分之一
                end_pos = round((len(seq) + 2 )/(frag_num_left+1))
                print("end_pos:", end_pos)
                frag, rest_seq, overhang, overhang_pair = split_with_pos(seq, end_pos, overhang_len)
                print(len(frag), len(rest_seq))
                #try_times = 0
                while not check_overhang_unique(overhang_list[i-1], overhang):
                    print(overhang)
                    if i == 1:
                        for j in range(end_pos, max_len_5 + 1):
                            print(j)
                            frag, rest_seq, overhang, overhang_pair = split_with_pos(seq, j, overhang_len)
                            if check_overhang_unique(overhang_list[i-1], overhang):
                                print(j, overhang)
                                print("pick end_pos from end_pos to max_len_5")
                                break
                    else:
                        for j in range(end_pos, max_inner_len + 1):
                            print(j)
                            frag, rest_seq, overhang, overhang_pair = split_with_pos(seq, j, overhang_len)
                            if check_overhang_unique(overhang_list[i-1], overhang):
                                print(j, overhang)
                                print("pick end_pos from end_pos to max_inner_len")
                                break
                    break
                while not check_overhang_unique(overhang_list[i-1], overhang):
                    print(overhang)
                    for j in range(end_pos, frag_min_len - 1, -1):
                        print(j)
                        frag, rest_seq, overhang, overhang_pair = split_with_pos(seq, j, overhang_len)
                        if check_overhang_unique(overhang_list[i-1], overhang):
                            print(j, overhang)
                            print("pick end_pos from frag_min_len to end_pos")
                            break
                    break  
                while not check_overhang_unique(overhang_list[i-1], overhang):
                    print(overhang_list[i-1], overhang)
                    print('Oligo %s can\'t find unique cohesive end.' %name)
                    break

                    '''end_pos = np.random.randint(frag_min_len, max_len_5)
                    try_times += 1
                    frag, rest_seq, overhang, overhang_pair = split_with_pos(seq, end_pos, overhang_len)
                    if try_times == 100:
                        print('Oligo %s can\'t find unique cohesive end.' %name)
                        print("Try again!")
                        break'''
                # if try_times == 100:
                #     break
                overhang_list[i-1].append(overhang)
                overhang_list[i-1].append(overhang_pair)
                seq = rest_seq
                if i == 1:                 
                    #添加adapters
                    frag_1 = add_adapter(frag, ecoli_primer5, oligoa_3[batch_big_id], bsai)
                    
                    split_result['a_'+name] = frag_1
                    
                elif frag_num_left == 1:            
                    #用于添加adapters
                    frag_2 = add_adapter(frag, oligob_5[batch_big_id], oligob_3[batch_big_id], bsmi)
                    frag_3 = add_adapter_last(rest_seq, oligoc_5[batch_big_id], ecoli_primer3, BsmI)
                    #将每一条序列赋予新的命名
                    split_result['b_'+name] = frag_2
                    split_result['c_'+name] = frag_3      
        if not check_overhang_unique(overhang_list[i-1], overhang):
            continue
        else:
            print("Work is done!\n")
            cut_successful_flag = True
    
    return split_result, overhang_list

def main(args):
    ####################################
    #input oligo length
    max_oligo_size = args.max_oligo_size
    min_len = args.min_oligo_size
    frag_num = args.frag_num
    ad_fname = args.adapter_fname
    adapter_idx = args.adapter_number - 1

    #input the design names and sequences
    with open(args.input_list,'r', encoding='utf-8-sig') as ip:
        lines = ip.readlines()

    seq_list = {}
    for line in lines:
        its = line.strip().split()
        if len(its) != 2:
            print("Skipping:")
            print(its)
            continue
        else:
            seq_list[its[0]] = its[1]

    print(f'file loaded! {len(seq_list)}')

    '''
    BEGGINING OF MAIN:
    '''
    # Load and process sequence files:
    # for this verison, specify which line of the adapters should be used
    # numbering starts from 1, not like python 0
    # Get adapters ready
    #oligoa_5, oligoa_3, oligob_5, oligob_3, oligoc_5, oligoc_3 = read_adapters(ad_fname, frag_num)#修改函數 返回列表
    adapter_list = read_adapters(ad_fname, frag_num)


    #split the sequences into oligos
    result_list, overhang_list = split_sequences_more_fragments(seq_list, frag_num, adapter_list, adapter_idx)
    
    #write the output file
    with open('%d_oligos_subpools_%d.tab'%(frag_num, args.adapter_number), 'w') as outputfile:
        for key, value in result_list.items():
            if len(value) > 300:
                print('Oligo %s is longer than 300bp!!!' % key)
            if len(value) < 250:
                print('Oligo %s is shorter than 250bp!!!' % key)
            outputfile.write(key + ', ' + value + '\n')


if __name__ == '__main__':
    ########## Option system: ###########
    argparser = argparse.ArgumentParser(description='Split genes in orthogonal pieces that can be used in multiplex assembly')
    argparser.add_argument('-input_list', type=str, help='Name of file containing: mygenename DNAsequence ')
    argparser.add_argument('-adapter_fname', type=str, default='./pool_adapters_short_list.txt',help='Name of file containing adapter sequences. Format: First line: column names, followed by lines: adapter_name fiveprime_5 fiveprime_3 threeprime_5 threeprime_3')
    argparser.add_argument('-adapter_number', type=int, help='What adapter to use? starting at 1')
    argparser.add_argument('-max_oligo_size', type=int, default=300, help='Absolute max length of orderable oligo')
    argparser.add_argument('-min_oligo_size', type=int, default=120, help='Absolute min length of orderable oligo')
    argparser.add_argument('-enzyme_num', type=int, default=1, choices=[1, 2], help='The number of how many enzymes.')
    argparser.add_argument('-frag_num', type=int, default=3, help='The number of how many fragments you want to split.')
    
    args = argparser.parse_args()
    main(args)

