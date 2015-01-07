#!/usr/bin/python
# Convert DNA sequence to its complementary sequence
# Input file should be in fasta format
# Replace "complement_seq" with "reverse_complement_seq" in 
# return statement if reverse complementary sequence is needed

base_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

def convert_seq(seq):
    complement_seq = ''
    for base in seq:
        complement_seq += base_dict[base.upper()]
    reverse_complement_seq = complement_seq[::-1]
    return complement_seq

with open('c:\seq_output.txt', 'w') as output_seq:
    with open('c:\seq.txt', 'r') as input_seq:
        for line in input_seq:
            if line.startswith('>'):
                print >> output_seq, line.strip()
            else:
                print >> output_seq, line.strip().upper()
                print >> output_seq, convert_seq(line.strip())
