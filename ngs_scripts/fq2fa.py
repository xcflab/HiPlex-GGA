from Bio import SeqIO
import sys

# 
input_file = sys.argv[1]
output_file = sys.argv[2]

#
with open(output_file, "w") as output_handle:
    SeqIO.convert(input_file, "fastq", output_handle, "fasta")
