#!/usr/bin/env python
# Aug 28, 2024
# Zhien Wu modified

import numpy as np
from tqdm import tqdm
import sys
import random
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from multiprocessing import Pool

#overhand_list
#The best overhangs are included in overhangs_good. The second good overhangs are included in bad_overhangs.
overhangs_good = ['TCAA', 'ATAA', 'TAGA', 'GTTA', 'AGGG', 'CCTA', 'AAGA', 'TCCA', 'AGAA', 'AAAG', 'ACAT', 'GGGA', 'GCAA', 'TAAA', 'TGAA', 'GACA', 'ACGA', 'GGTA', 'ATCC', 'ATTG', 'CTAC', 'AAGC', 'CATC', 'ACTC', 'CACA', 'CTAA', 'GAAA', 'AGAC', 'AGCA', 'CGTA', 'ACCG', 'CAGA', 'AACT', 'AATA', 'GCAC', 'CCAG', 'CAAG', 'AAAT', 'ATCA', 'CAGG', 'CATA', 'GGAA', 'AGGA', 'ACGC', 'ATAC', 'CTCA', 'GCCA', 'CCGA', 'ACAG', 'AATC', 'CAGC', 'AAAA', 'AGTG', 'CGCA', 'AACG', 'GAGA', 'ACTA', 'TACA', 'ATGA', 'CGAC', 'CGAA', 'AGCC']
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

class Frag(object):
    """
    This object is used to store the information of the fragment.
    """
    def __init__(self, frag_num, dna_seq, aa_seq, max_len_5, max_inner_len, max_len_3, min_len): 
        self.frag_num = frag_num
        self.dna_seq = dna_seq
        self.max_len_5 = max_len_5
        self.max_inner_len = max_inner_len
        self.max_len_3 = max_len_3
        self.min_len = min_len
        self.aa_seq = aa_seq
        self.init_dna_seq = dna_seq
        self.overhang_len = 4
        self.end_pos = 0
        self.cut_order = 0
        self.frag_min_len = None
        
    def count_frag_lenth(self, cut_order, min_len=150):
        #This function is used to calculate the minimum length of the fragment.
        seq_len = len(self.dna_seq)
        frag_num_left = self.frag_num - self.cut_order
        frag_min_len = seq_len - (self.max_inner_len - 4) * ( frag_num_left - 1 ) - ( self.max_len_3 - 4)

        if frag_min_len < min_len:
            frag_min_len = min_len

        if self.cut_order == 1:
            assert frag_min_len <= self.max_len_5
        else:
            assert frag_min_len <= self.max_inner_len

        # if frag_num_left == 1:  #haven't used
        #     last_len = seq_len - self.max_inner_len
        #     if last_len < self.min_len:
        #         self.max_inner_len = seq_len - self.min_len * ( frag_num_left - 1 )
                
        return frag_min_len
    
    def update_dna_seq(self):
        new_dna_seq_begin = len(self.init_dna_seq) - len(self.dna_seq)
        self.dna_seq = self.init_dna_seq[new_dna_seq_begin:]
    
#function:
def sort_by_length(input_dict):
    """
    Sort the dictionary by the length of the values.
    """
    sorted_dict = {k: v for k, v in sorted(input_dict.items(), key=lambda item: len(item[1]),reverse=True)}
    return sorted_dict

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

def check_overhang_unique(end_list, to_check_end):
    """
    Check whether the overhang is unique
    """
    flag = True
    if to_check_end in end_list:
        flag = False
    return flag

def split_with_pos(sequence, pos, overhang_len=4):
    """
    Split the sequence with the position
    """
    end_len_1 = overhang_len
    seq_1 = sequence[0:pos]
    seq_2 = sequence[pos-end_len_1:]
    overhang = Seq(sequence[pos-end_len_1:pos])
    overhang_pair = overhang.reverse_complement()
    return seq_1, seq_2, overhang, overhang_pair

def get_frame_end(startpoint):
    """
    Returns the frame end given the that the sequence has length "startpoint"
    """
    frameend = startpoint

    if frameend % 3 == 0 :
        frameend -= 1
    elif frameend % 3 == 1:
        frameend += 1

    return frameend

def replace_codons(s,frame_end,overlap_len,codons):
    """
    Returns a new sequence encoding for the same peptide, but with different codons in the overlap region.
    """
    #starting at the first position for new codon hence +1
    randomized_positions = [i for i in range((frame_end+1)//3, (frame_end+1+ overlap_len)//3 )] #LA / to // to make integer
    random.shuffle(randomized_positions)

    #change 6 codons at a time
    poss_to_change = min(6,len(randomized_positions)-1)
    changed = randomized_positions[:poss_to_change]
    a = sorted(changed)

    sh = ""
    for p in range(frame_end+1,frame_end+1+overlap_len, 3):
        #if is_part(a,p/3):
        if p/3 in a:
            #Make a new codon
            codon = ''
            for i in range(3):
                codon += s[p + i]
            #grab a codon at random
            aa = str(Seq(codon).translate())
            r = random.randint(1,len(codons[aa]))
            sh += codons[aa][r-1]
        else:
            # Return original codon
            for i in range(3):
                sh += s[p + i]
    return sh

def find_overhang_from_appropriate_range(gene, overhang_list, trial = -1):
    """
    This function is used to find the unique overhang from the appropriate range.
    """
    reasonable_tag = True
    
    max_len = gene.max_len_5 if gene.cut_order == 1 else gene.max_inner_len
    if max_len < gene.end_pos:
        reasonable_tag = False
        return reasonable_tag
    else:
        for j in range(gene.end_pos, max_len + 1):
            frag, rest_seq, overhang, overhang_pair = split_with_pos(gene.dna_seq, j, gene.overhang_len)
            if overhang in overhangs_good and check_overhang_unique(overhang_list, overhang):
                break
            if trial > 499:
                if overhang not in bad_overhangs and check_overhang_unique(overhang_list, overhang):
                    break
        
        #If there is no unique overhang from end_pos to max_len, return the last overhang
        return frag, rest_seq, overhang, overhang_pair
    
def find_overhang(gene, overhang_list, codons, overhang_len = 4):
    """
    This function is used to find the unique overhang for the fragment.
    """
    
    split_successful_tag = False
    # If the 5' end fragment is shorter than the inner fragments, and this is the first cutting operation,
    # and end_pos is greater than or equal to the maximum length of the 5' end fragment,
    # then adjust end_pos to the minimum fragment length to ensure enough space for the 5' end fragment.
    if gene.max_len_5 < gene.max_inner_len and gene.cut_order == 1 and gene.end_pos >= gene.max_len_5:
        gene.end_pos = gene.frag_min_len
        
    # In order to ensure the rest fragments have enough space, adjust end_pos to the minimum fragment length.
    if gene.end_pos < gene.frag_min_len:
        gene.end_pos = gene.frag_min_len
        
    # Find the unique overhang
    while not split_successful_tag:
        result = find_overhang_from_appropriate_range(gene, overhang_list)
        if isinstance(result, bool):
            break
        else:
            frag, rest_seq, overhang, overhang_pair = result
        
        #Check whether the overhang is unique and good enough, if not, redesign the dna sequence
        if overhang in overhangs_good and check_overhang_unique(overhang_list, overhang):
            split_successful_tag = True
        else:
            for trial in range(1000):
                #Check whether the overhang is acceptable
                if overhang in overhangs_good and check_overhang_unique(overhang_list, overhang):
                    split_successful_tag = True
                    break
                if trial > 499:
                    if overhang not in bad_overhangs and check_overhang_unique(overhang_list, overhang):
                        split_successful_tag = True
                        break 
                #If the overhang is not accepted, redesign the dna sequence
                split_posi = len(gene.init_dna_seq) - len(gene.dna_seq) + gene.end_pos # the position of redesign dna sequence begins
                fiveprime = gene.init_dna_seq[:split_posi] # the 5' end of the dna sequence
                frame_end = get_frame_end(len(fiveprime)) # the frame end of the dna sequence
                change_len = gene.max_len_5 - gene.end_pos if gene.cut_order == 1 else gene.max_inner_len - gene.end_pos # the length of the dna sequence that needs to be redesigned
                newseq = replace_codons(gene.init_dna_seq,frame_end,change_len,codons) # the redesigned dna sequence
                gene.init_dna_seq = gene.init_dna_seq[:frame_end+1] + newseq + gene.init_dna_seq[frame_end+1+len(newseq):] # the new dna sequence
                gene.update_dna_seq() # update the dna sequence
                #Find the unique overhang
                result = find_overhang_from_appropriate_range(gene, overhang_list, trial)
                if isinstance(result, bool):
                    break
                else:
                    frag, rest_seq, overhang, overhang_pair = result
            #After 1000 trials, if the overhang is still not unique, break the loop
            if not check_overhang_unique(overhang_list, overhang):
                break
            
    if split_successful_tag:
        return frag, rest_seq, overhang, overhang_pair
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
    """
    Add the spool barcode to the fragments.
    """
    
    alphabet = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z"]
    frag_ad_list = {}
    
    for i in range(frag_num):
        name = alphabet[i]+ "_"+ seq_name
        if frag_seq[i] == "":
            frag_ad_list[name] = "not_split"
        else:
            if i == 0:
                # the first fragment
                if seq_barcode != '':
                    seq_barcode = NcoI + seq_barcode[0]
                
                before_compelment = adapter_F + spool_barcode_list[2 * i][subp_barc_idx] + seq_barcode + NdeI + frag_seq[i] + bsai
                after_compelment = spool_barcode_list[2*i + 1][subp_barc_idx] + adapter_R 
                
                frag_ad = before_compelment + after_compelment
                if len(frag_ad) < 211 :
                    frag_ad = before_compelment + 'atgagccatattcaacgggaaacgtcttgctgcgattaaattccaacatggatgctgatttatatgggtatataat' + bsai + after_compelment
                elif len(frag_ad) < 251 :
                    frag_ad = before_compelment + 'atgagccatattcaacgggaaacgtcttgctgtaat' + bsai + after_compelment
                
                frag_ad_list[name] = frag_ad
                
            elif i == frag_num - 1:
                # the last fragment
                if i % 2 == 0:
                    enzyme = BsmBI
                else:
                    enzyme = BsaI
                if seq_barcode != '':
                    seq_barcode = seq_barcode[1] + BamHI 
                    
                before_compelment = adapter_F + spool_barcode_list[2 * i][subp_barc_idx] + enzyme
                after_compelment = frag_seq[i] + XhoI + seq_barcode + spool_barcode_list[2*i + 1][subp_barc_idx] + adapter_R
                
                frag_ad = before_compelment + after_compelment
                if len(frag_ad) < 211:
                    frag_ad = before_compelment + 'atgagccatattcaacgggaaacgtcttgctgcgattaaattccaacatggatgctgatttatatgggtatataat' + enzyme + after_compelment
                elif len(frag_ad) <251:
                    frag_ad = before_compelment + 'atgagccatattcaacgggaaacgtcttgctgtaat' + enzyme + after_compelment
                
                frag_ad_list[name] = frag_ad
            
            else:
                # the inner fragments
                if i % 2 == 0:
                    enzyme_5 = BsmBI
                    enzyme_3 = bsai
                else:
                    enzyme_5 = BsaI
                    enzyme_3 = bsmbi
                
                before_compelment = adapter_F + spool_barcode_list[2 * i][subp_barc_idx] + enzyme_5 + frag_seq[i] + enzyme_3
                after_compelment = spool_barcode_list[2*i + 1][subp_barc_idx] + adapter_R
                
                frag_ad = before_compelment + after_compelment
                if len(frag_ad) < 211:
                    frag_ad = before_compelment + 'atgagccatattcaacgggaaacgtcttgctgcgattaaattccaacatggatgctgatttatatgggtatataat' + enzyme_3 + after_compelment
                elif len(frag_ad) < 251 :
                    frag_ad = before_compelment + 'atgagccatattcaacgggaaacgtcttgctgtaat' + enzyme_3 + after_compelment
                
                frag_ad_list[name] = frag_ad
 
    return frag_ad_list

def split_sequences(designs, 
                    protein_seqs,
                    codons,
                    frag_num, 
                    spool_barcode_list, 
                    subp_barc_idx, 
                    adapter_F, 
                    adapter_R, 
                    seq_barcode_list, 
                    max_oligo_size=300,
                    min_oligo_len=250):
    """
    This function is used to split the sequences into oligos.
    """
    #Check if the size of subpool is larger than 100
    assert len(designs) < 100
    
    #subpool barcode length, all the same
    subp_barc_len_5 = len( spool_barcode_list[0][subp_barc_idx] ) #  the length of 5' subpool barcode
    subp_barc_len_3 = len( spool_barcode_list[-1][subp_barc_idx] ) #  the length of 3' subpool barcode
    subp_barc_in_len = len( spool_barcode_list[1][subp_barc_idx] ) #  the length of inner subpool barcode
    #adapter length
    adapter_F_len = len( adapter_F ) #  the length of 5' adapter
    adapter_R_len = len( adapter_R ) #  the length of 3' adapter
    #sequence barcode length
    if seq_barcode_list == None:
        seq_barc_len = 0
    else:
        seq_barc_len = len(seq_barcode_list[0][0]) + 6 #  the length of sequence barcode

    #the length of the longest A fragment
    max_len_5 = max_oligo_size - adapter_F_len - subp_barc_len_5 - 6 - seq_barc_len - (subp_barc_in_len + 7) - adapter_R_len #  add the basi or not?
    #the length of the longest inner fragment
    max_inner_len = max_oligo_size - adapter_F_len - (subp_barc_in_len + 7) * 2 - adapter_R_len
    #the length of the longest last fragment
    max_len_3 = max_oligo_size - adapter_F_len - (subp_barc_in_len + 7) - 6 - seq_barc_len - subp_barc_len_3 - adapter_R_len
    
    #generate overhang list
    overhang_list = []
    for frag_order in range(1, frag_num):
        frag_order = []
        overhang_list.append(frag_order)

    split_result = {}
    
    #split the sequences
    seq_num = 0
    
    designs = sort_by_length(designs)
    for name, init_dna_seq in tqdm(designs.items()):
        # Sequence initialization
        aa = protein_seqs[name]
        gene = Frag(frag_num=frag_num, dna_seq=init_dna_seq, max_len_5=max_len_5, max_inner_len=max_inner_len, max_len_3=max_len_3, min_len=min_oligo_len, aa_seq=aa)

        cut_successful_flag = False
        not_unique_overhang = False
        
        try_times = 0
        while not cut_successful_flag:
            frag_seq = []
            gene.dna_seq = gene.init_dna_seq
            for i in range(1, frag_num):
                gene.cut_order = i
                frag_num_left = frag_num - i
                gene.frag_min_len = gene.count_frag_lenth(i)
                gene.end_pos = round((len(gene.dna_seq) + frag_num_left * 4 )/(frag_num_left+1))#End_pos is the average length of each fragment.

                #find unique overhang
                result = find_overhang(gene, overhang_list[i-1], codons, overhang_len = 4)
                
                if isinstance(result, bool):
                    #can't find unique overhang, need to re-split
                    try_times += 1
                    break
                else:
                    frag, rest_seq, overhang, overhang_pair = result
                    gene.dna_seq = rest_seq
                    #Put frag into a list 
                    frag_seq.append(frag)
                    if i == frag_num - 1:
                        frag_seq.append(rest_seq)
                    #Add overhang to the list
                    overhang_list[i-1].append(overhang)
                    overhang_list[i-1].append(overhang_pair) 

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
            # optional: add sequence barcode
            if seq_barcode_list == None:
                seq_barcode = ''
            else:
                seq_barcode = seq_barcode_list[seq_num]
           
           # Add adapters and spool barcodes 
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
    if len(split_result) < len(designs) or not_suitable_length_flag:
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
            
    designs = {}
    protein_seqs = {}
    for line in lines:
        its = line.strip().split()
        if len(its) != 3:
            print("Skipping:")
            print(its)
            continue
        else:
            #input format name\tAA_seq\tNucleotide_seq
            protein_seqs[its[0]] = its[1]
            designs[its[0]] = its[2]
            
    print(f'file loaded! {len(protein_seqs)}')
    
    #Load codons
    codonfile = open(args.codontable_fname, 'r').readlines()
    codons = {}
    for line in codonfile:
        its = line.strip().split()
        try:
            if not codons[its[0]]:
                codons[its[0]].append(its[1])

            else:
                codons[its[0]].append(its[1])
        except KeyError:
            bases = []
            codons[its[0]] = bases
            codons[its[0]].append(its[1])

    print("your condons\n", codons)
    
    '''
    BEGGINING OF MAIN:
    '''
    # Load and process sequence files:
    # for this verison, specify which line of the adapters should be used
    # numbering starts from 1, not like python 0
    
    # Get adapters ready
    spool_barcode_list = read_subpool_barcode(spool_barc_fname, frag_num)
    
    # Get sequence barcodes ready
    if seq_barc_fname == None:
        seq_barcode_list = None
    else:
        seq_barcode_list = read_sequence_barcode(seq_barc_fname, len(designs))
    
    #split the sequences into oligos
    result = split_sequences(designs, 
                             protein_seqs,
                             codons, 
                             frag_num, 
                             spool_barcode_list, 
                             subp_barc_idx, 
                             adapter_F, adapter_R, 
                             seq_barcode_list, 
                             max_oligo_size = args.max_oligo_length,
                             min_oligo_len = args.min_oligo_length)
    
    while isinstance(result, bool):
        designs = shuffle_dict(designs)
        result = split_sequences(designs, 
                                protein_seqs,
                                codons, 
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
    argparser.add_argument('--input_list', type=str, required=True, help='Name of file containing: mygenename DNAsequence ')
    argparser.add_argument('--subpool_barcode_fname', type=str, default='./pool_subpools_barcode.txt',help='Name of file containing adapter sequences. Format: First line: column names, followed by lines: adapter_name fiveprime_5 fiveprime_3 threeprime_5 threeprime_3')
    argparser.add_argument('--adapter_f', type=str, default='FFFFFFFFFFFFFFFFFFFF',help='Forward adapter sequences. Default: ')
    argparser.add_argument('--adapter_r', type=str, default='RRRRRRRRRRRRRRRRRRRR',help='Reverse adapter sequences. Default: ')
    argparser.add_argument('--sequence_barcode_fname', type=str, help='Name of file containing sequence barcode sequences. Format: First line: column names, followed by lines: fiveprime_5 threeprime_3')
    argparser.add_argument('--subp_barc_index', type=int, required=True,  help='What subpool barcode to use? starting at 1')
    argparser.add_argument('--max_oligo_length', type=int, default=300, help='Absolute max length of orderable oligo')
    argparser.add_argument('--min_oligo_length', type=int, default=251, help='Absolute min length of orderable oligo')
    argparser.add_argument('--frag_num', type=int, default=3, help='The number of how many fragments you want to split.')
    argparser.add_argument('-codontable_fname', type=str,default='./codontable.tab',help='Codon table to use')
    argparser.add_argument('-nproc', type=int,default=1,help='Number of processors to use. Must be in the same node (DIGs: -cX -N1 AFAIK)')
    
    args = argparser.parse_args()
    main(args)

