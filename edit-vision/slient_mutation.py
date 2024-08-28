#!/usr/bin/env python
# Mar 12, 2024
# Zhien Wu modified
import argparse
from tqdm import tqdm

'''
This script can help the user mutate the enzyme site. 
'''


    
#functions
    
def change_typeii_5(site, seq, name):
    enzymes = {1: 'C', 2: 'C', 0: 'T'}  # 余数1对应的替换字符为C，余数2对应的替换字符为C，余数0对应的替换字符为T
    offset = {1: 2, 2: 4, 0: 3}  
    
    print(f"{name} has a BsaI enzyme site in {site + 1}") 
    site_seq = seq[site:site+6]  
    remainder = (site + 1) % 3  
    
    if remainder in enzymes:  
        replace_char = enzymes[remainder]  
        replace_pos = site + offset[remainder]  
        seq = seq[:replace_pos] + replace_char + seq[replace_pos + 1:]  
        print(f"Change {site_seq} into {seq[site:site+6]}")  
        
    assert site_seq != seq[site:site+6]  
    return seq  

def change_typeii_3(site, seq, name):
    enzymes = {1: 'A', 2: 'T', 0: 'C'}  # 余数1对应的替换字符为A，余数2对应的替换字符为T，余数0对应的替换字符为C
    offset = {1: 2, 2: 4, 0: 1} 
    
    print(f"{name} has a BsaI enzyme site in {site + 1}") 
    site_seq = seq[site:site+6] 
    remainder = (site + 1) % 3 
    
    if remainder in enzymes:  
        replace_char = enzymes[remainder]  
        replace_pos = site + offset[remainder]  
        seq = seq[:replace_pos] + replace_char + seq[replace_pos + 1:] 
        print(f"Change {site_seq} into {seq[site:site+6]}") 
        
    assert site_seq != seq[site:site+6] 
    return seq  

def change_typei(site, seq, name, pattern):
    if pattern == 'CATATG':  
        enzyme = {1: 'A', 2: 'T', 0: 'C'}
        offset = {1: 2, 2: 4, 0: 1}
        print(f"{name} has a BsaI enzyme site in {site + 1}") 
        site_seq = seq[site:site+6]  
        remainder = (site + 1) % 3  
    
        if remainder in enzymes:  
            replace_char = enzymes[remainder]  
            replace_pos = site + offset[remainder]  
            seq = seq[:replace_pos] + replace_char + seq[replace_pos + 1:]  
            print(f"Change {site_seq} into {seq[site:site+6]}")  
    elif pattern == 'CTCGAG':  
        replace_char = 'C'  
        replace_pos = site + 4
    elif pattern == 'AAGCTT':  
        replace_char = 'C'  
        replace_pos = site + 4
    elif pattern == 'GAATTC':  
        replace_char = 'C'  
        replace_pos = site + 4
    elif pattern == 'CCATGG':
        replace_char = 'C'
        replace_pos = site + 4
    elif pattern == 'GGATCC':
        replace_char = 'C'
        replace_pos = site + 4
    
    return seq
    
def change_bsai5(seq, name):
    seq5_bsai = 'GGTCTC'
    while seq5_bsai in seq:  
        site = seq.find(seq5_bsai)  
        seq = change_typeii_5(site, seq, name)  
    return seq

def change_bsai3(seq, name):
    seq3_bsai = 'GAGACC'
    while seq3_bsai in seq:  
        site = seq.find(seq3_bsai)  
        seq = change_typeii_3(site, seq, name)  
    return seq

def change_BsmBI5(seq, name):
    seq5_bsmbi = 'CGTCTC'
    while seq5_bsmbi in seq:  
        site = seq.find(seq5_bsmbi)  
        seq = change_typeii_5(site, seq, name)  
    return seq

def change_BsmBI3(seq, name):
    seq3_bsmbi = 'GAGACG'
    while seq3_bsmbi in seq:  
        site = seq.find(seq3_bsmbi)  
        seq = change_typeii_3(site, seq, name) 
    return seq

def change_the_rest(seq, name):
    sites = {
        'NdeI': 'CATATG',
        'XhoI': 'CTCGAG',
        'HindIII': 'AAGCTT',
        'EcoRI': 'GAATTC',
        'NcoI': 'CCATGG',
        'BamHI': 'GGATCC'
    }
    
    for enzyme, pattern in sites.items():
        while pattern in seq:
            site = seq.find(pattern)
            seq = change_typei(site, seq, name, pattern)
    return seq

def translate(seq):
    """
    Translate a string containing a nucleotide sequence into a string containing the corresponding sequence of amino acids .
    Nucleotides are translated in triplets using the table dictionary; each amino acid 4 is encoded with a string of length 1.
    """
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  #this table dictionary is pre-created
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',

        'ata':'i', 'atc':'i', 'att':'i', 'atg':'m',
        'aca':'t', 'acc':'t', 'acg':'t', 'act':'t',
        'aac':'n', 'aat':'n', 'aaa':'k', 'aag':'k',
        'agc':'s', 'agt':'s', 'aga':'r', 'agg':'r',                  #this table dictionary is pre-created
        'cta':'l', 'ctc':'l', 'ctg':'l', 'ctt':'l',
        'cca':'p', 'ccc':'p', 'ccg':'p', 'cct':'p',
        'cac':'h', 'cat':'h', 'caa':'q', 'cag':'q',
        'cga':'r', 'cgc':'r', 'cgg':'r', 'cgt':'r',
        'gta':'v', 'gtc':'v', 'gtg':'v', 'gtt':'v',
        'gca':'a', 'gcc':'a', 'gcg':'a', 'gct':'a',
        'gac':'d', 'gat':'d', 'gaa':'e', 'gag':'e',
        'gga':'g', 'ggc':'g', 'ggg':'g', 'ggt':'g',
        'tca':'s', 'tcc':'s', 'tcg':'s', 'tct':'s',
        'ttc':'f', 'ttt':'f', 'tta':'l', 'ttg':'l',
        'tac':'y', 'tat':'y', 'taa':'_', 'tag':'_',
        'tgc':'c', 'tgt':'c', 'tga':'_', 'tgg':'w',
    }
    protein=""
    if (len(seq)-1)%3==0:
        for i in range(0,len(seq)-1,3):
            codon=seq[i:i+3]
            protein+=table[codon]
    if (len(seq)-2)%3==0:
        for i in range(0,len(seq)-2,3):
            codon=seq[i:i+3]
            protein+=table[codon]
    if len(seq)%3==0:
        for i in range(0,len(seq),3):
            codon=seq[i:i+3]
            protein+=table[codon]
    return protein

def read_seq(inputfile):
    with open(inputfile,"r") as f:
        seq=f.read()
    seq=seq.replace("\n","")
    seq=seq.replace("\r","")
    return seq

def main(args):
    '''
    BEGGINING OF MAIN:
    '''
    #input the sequences and names
    with open(args.input_list, 'r') as ip:
        lines = ip.readlines()

    seq_list = {}
    for line in lines:
        its = line.strip().split()
        if len(its) != 2:
            seq_list[its[0]] = its[2]
        else:
            seq_list[its[0]] = its[1]

    print(f'file loaded! {len(seq_list)}')

    changed_seq_list = {}

    enzyme_functions = {
    'BsaI': (change_bsai5, change_bsai3),
    'BsmBI': (change_BsmBI5, change_BsmBI3),
    'Both': (change_bsai5, change_bsai3, change_BsmBI5, change_BsmBI3),
    'All': (change_bsai5, change_bsai3, change_BsmBI5, change_BsmBI3, change_the_rest)
    }

    enzyme_actions = enzyme_functions.get(args.enzyme_type)

    if enzyme_actions:
        # 执行选择的修改操作
        print(f'Change the {args.enzyme_type} site(s).')
        for name, seq in tqdm(seq_list.items()):
            beforechange = translate(seq)
            for action in enzyme_actions:
                seq = action(seq, name)
            afterchange = translate(seq)
            assert beforechange == afterchange
            changed_seq_list[name] = seq
    else:
        print("Invalid enzyme type.")

    #write the output file
    with open('change_%s.txt'%args.enzyme_type, 'w') as outputfile:
        for key, value in changed_seq_list.items():
            outputfile.write(key + ' ' + value + '\n')
            
if __name__ == '__main__':
    ########## Option system: ###########
    argparser = argparse.ArgumentParser(description='Slient mutation of enzyme site')
    argparser.add_argument('-input_list', type=str, help='Name of file containing: mygenename (AAsequence) DNAsequence ')
    argparser.add_argument('-enzyme_type', type=str, default='All', choices=['BsaI', 'BsmBI', 'Both', 'All'], help='enzyme type: BsaI, BsmBI, or Both; NdeI, XhoI, HindIII, EcoRI or all')
    
    args = argparser.parse_args()
    main(args)

