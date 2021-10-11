#!/bin/usr/python 3

__author__ = "Jurrien de Jong"
__date__ = "8/10/2021"
__version__ = "V1.0"

from Bio import SeqIO
import sys
import argparse
from Bio import AlignIO

test = 5

def align(file_name):
    data = []
    if file_name.endswith(".msf"):
        alignment = AlignIO.read(open(file_name), "msf")
        for record in alignment:
            data.append(record.seq)
        return data
    else:
        raise FileNotFoundError


def dna_to_protein(sequence):
    protein_seq = ""
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
    if len(sequence) % 3 == 0:
        for i in range(0, len(sequence), 3):
            codon = sequence[i:i + 3]
            protein_seq += table[codon]
    return protein_seq


def main(arguments):
    data = align(arguments.in_file)
    protein_seq = dna_to_protein(arguments.Sequence)
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process a MSA.')
    parser.add_argument('-i', type=str, dest='in_file', help='Please hand over a .msf file')
    parser.add_argument('-seq', type=int, dest='Sequence', help='Please hand over a DNA or protein sequence.')
    parser.add_argument('-pos', type=int, dest='SNP_Pos', help='Please hand over a SNP position in the sequence.')
    parser.add_argument('-snp', type=str, dest='SNP', help='Please hand over a valid SNP.')
    arguments = parser.parse_args()
    sys.exit(main(arguments))
