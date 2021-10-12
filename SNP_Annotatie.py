#!/bin/usr/python 3

__author__ = "Jurrien de Jong"
__date__ = "8/10/2021"
__version__ = "V1.0"

import re
import sys
import argparse
from Bio import AlignIO
from Bio import Align


def read_file(file_name):
    data = []
    if file_name.endswith(".msf"):
        alignment = AlignIO.read(open(file_name), "msf")
        for record in alignment:
            data.append(record.seq.strip("\n"))
        return data
    else:
        raise Exception("An incorrect filetype is given."
                        " Please hand over a .msf file")


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
    if not re.match("[BDEFHIJKLMNOPQRSUVWXYZ]", sequence):
        for i in range(0, len(sequence), 3):
            codon = sequence[i:i + 3]
            protein_seq += table[codon]
        print("Using protein sequence : {}".format(protein_seq))
        return protein_seq
    else:
        return sequence


def implement_snp(seq, pos, snp):
    if len(snp) > 1:
        raise Exception("SNP should be one char long")
    if pos in range(0, len(seq)):
        sequence = seq[:pos] + snp + seq[pos + 1:]
    else:
        raise Exception("This index does not exist"
                        " because the sequence is "
                        "length: {}".format(len(seq) - 1))
    return sequence


def get_align_score(msa, snp_seq):
    scores = []
    aligner = Align.PairwiseAligner()
    for x in range(0, len(snp_seq)):
        for seq in msa:
            alignments = aligner.align(seq[x], snp_seq[x])
            scores.append(alignments.score)

    added_scores = [sum(scores[x:x+20]) for x in range(0, len(snp_seq) * len(msa), 20)]
    return added_scores


def write_results(scores, pos):
    message = "Your score is not valid!"
    if scores[pos] <= 2:
        message = "So this is a pretty bad SNP, it has a lot of consequences!\n" \
                  "The AA might be pretty preserved!"
    if 2 < scores[pos] <= 7:
        message = "So this SNP has some bad affects, but is not terrible."
    if 7 < scores[pos] <= 13:
        message = "So the SNP has very few bad effect."
    if scores[pos] > 13:
        message = "So this SNP has a neutral effect! :)"
    print("Severity of SNP at pos: {}, has score: {}.\n{}"
          .format(pos, scores[pos], message))


def main(arguments):
    print("\n\n\n- Starting Program -\n")
    msa = read_file(arguments.in_file)
    protein_seq = dna_to_protein(arguments.Sequence)
    snp_seq = implement_snp(protein_seq,
                            arguments.SNP_Pos,
                            arguments.SNP)

    scores = get_align_score(msa, snp_seq)
    write_results(scores, arguments.SNP_Pos)
    print("\n- End of Program -")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process a MSA.')
    parser.add_argument('-i', type=str, dest='in_file', help='Please hand over a .msf file')
    parser.add_argument('-seq', type=str, dest='Sequence', help='Please hand over a DNA or protein sequence.')
    parser.add_argument('-pos', type=int, dest='SNP_Pos', help='Please hand over a SNP position in the sequence.')
    parser.add_argument('-snp', type=str, dest='SNP', help='Please hand over a valid SNP.')
    arguments = parser.parse_args()
    sys.exit(main(arguments))
