#!/usr/bin/env python
# Mar 28, 2024
# Zhien Wu modified
import argparse
from tqdm import tqdm

'''
This script can help the user mutate the enzyme site. 
'''


#----Functions----#
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

def find_new_codon(aa, old_codon):
    '''
    Find the new codon based on the amino acid and the old codon.
    '''
    #Modified accroding to the codon usage of E. coli
    #delete 'L: CTA', 'I: ATA', 'P: CCC', 'K: AAG', 'R: CGA CGG AGA AGG', 'G: GGA'
    table = {
        'I': ['ATC', 'ATT'],
        'M': ['ATG'],
        'T': ['ACA', 'ACC', 'ACG', 'ACT'],
        'N': ['AAC', 'AAT'],
        'K': ['AAA'],
        'S': ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'],
        'R': ['CGC', 'CGT'],
        'L': ['CTC', 'CTG', 'CTT', 'TTA', 'TTG'],
        'P': ['CCA', 'CCG', 'CCT'],
        'H': ['CAC', 'CAT'],
        'Q': ['CAA', 'CAG'],
        'V': ['GTA', 'GTC', 'GTG', 'GTT'],
        'A': ['GCA', 'GCC', 'GCG', 'GCT'],
        'D': ['GAC', 'GAT'],
        'E': ['GAA', 'GAG'],
        'G': ['GGC', 'GGG', 'GGT'],
        'F': ['TTC', 'TTT'],
        'Y': ['TAC', 'TAT'],
        '_': ['TAA', 'TAG', 'TGA'],
        'W': ['TGG'],
    }
    
    new_codon = table.get(aa)
    if new_codon:
        if old_codon in new_codon:
            new_codon.remove(old_codon)
            return new_codon[0]
        else:
            return new_codon[0]
    return None

def change_the_site(seq, name, enzyme_types_seq):
    sites = {
        'NdeI': 'CATATG',
        'XhoI': 'CTCGAG',
        'HindIII': 'AAGCTT',
        'EcoRI': 'GAATTC',
        'NcoI': 'CCATGG',
        'BamHI': 'GGATCC',
        'BsaI5': 'GGTCTC',
        'BsmBI5': 'CGTCTC',
        'BsaI3': 'GAGACC',
        'BsmBI3': 'GAGACG',
    }
    
    for enzyme, pattern in sites.items():
        while pattern in seq:
            site = seq.find(pattern)
            
            remainder = site % 3  #get the site of the frame
            if remainder == 0:
                old_codon = seq[site:site + 3]
            elif remainder == 1:
                if enzyme != 'NcoI':
                    old_codon = seq[site + 2:site + 5]
                else:
                    old_codon = seq[site - 1:site + 2]
            else:
                old_codon = seq[site + 1:site + 4]
                
            aa = translate(old_codon)
            new_codon = find_new_codon(aa, old_codon)
            
            if remainder == 0:
                seq = seq[:site] + new_codon + seq[site + 3:]
            elif remainder == 1:
                if enzyme != 'NcoI':
                   seq = seq[:site + 2] + new_codon + seq[site + 5:]
                else:
                    seq = seq[:site - 1] + new_codon + seq[site + 2:]
            else:
                seq = seq[:site + 1] + new_codon + seq[site + 4:]
                

    return seq

#----Main----#
def main(args):
    '''
    BEGGINING OF MAIN:  
    '''
    #input the sequences and names
    with open(args.input_list, 'r') as ip:
        lines = ip.readlines()

    filename = args.input_list.split('.')[0]
    
    seq_list = {}
    for line in lines:
        its = line.strip().split()
        if len(its) != 2:
            seq_list[its[0]] = its[2]
        else:
            seq_list[its[0]] = its[1]

    print(f'file loaded! {len(seq_list)}')

    changed_seq_list = {}

    enzyme_combinations = {
    'BsaI': ('BsaI5', 'BsaI3'),
    'BsmBI': ('BsmBI5', 'BsmBI3'),
    'Both': ('BsaI5', 'BsaI3', 'BsmBI5', 'BsmBI3'),
    'All': ('BsaI5', 'BsaI3', 'BsmBI5', 'BsmBI3', 'NdeI', 'XhoI', 'HindIII', 'EcoRI', 'NcoI', 'BamHI')
    }

    enzyme_types_seq = enzyme_combinations.get(args.enzyme_type)
    
    if enzyme_types_seq:
        for name, seq in tqdm(seq_list.items()):
            beforechange = translate(seq)
            seq = change_the_site(seq, name, enzyme_types_seq)
            afterchange = translate(seq)
            assert beforechange == afterchange
            changed_seq_list[name] = seq
        print("All enzyme sites have been changed.")
    else:
        print("Invalid enzyme type.")

    #write the output file
    with open('%s_change_%s.txt'%(filename, args.enzyme_type), 'w') as outputfile:
        for key, value in changed_seq_list.items():
            outputfile.write(key + ' ' + value + '\n')
            
if __name__ == '__main__':
    ########## Option system: ###########
    argparser = argparse.ArgumentParser(description='Slient mutation of enzyme site')
    argparser.add_argument('-input_list', type=str, help='Name of file containing: mygenename (AAsequence) DNAsequence ')
    argparser.add_argument('-enzyme_type', type=str, default='All', choices=['BsaI', 'BsmBI', 'Both', 'All'], help='enzyme type: BsaI, BsmBI, or Both; NdeI, XhoI, HindIII, EcoRI or all')
    
    args = argparser.parse_args()
    main(args)

