#!/bin/usr/python 3

"""
This program can accept a msa, sequence, SNP and it's position,
and hand out a message to the user if the SNP is 'bad'
"""

__author__ = "Jurrien de Jong"
__date__ = "8/10/2021"
__version__ = "V1.0"

import sys
import argparse
import math

from Bio import AlignIO
from Bio import Align


class SnpAnnotation:

    """
    This class can accept a msa, sequence, SNP and it's position,
    and hand out a message to the user if the SNP is 'bad'.
    It will create SNP annotation objects.
    """

    def __init__(self, file, pos, snp):
        self.file = file
        self.pos = pos
        self.snp = snp

    def read_msa(self):
        """

        Function which reads a .msf file

        :return alignment, column_len:
        """
        if self.file.endswith(".msf"):
            alignment = AlignIO.read(open(self.file), "msf")
        else:
            return Exception("\nPlease hand over a file with a .msf format!")

        # Get the length of the column
        column_len = len(alignment[:, 1])

        return alignment, column_len

    def process_sequence(self, seq_file):
        """
        Function which implements a SNP in a given seq,
        at a given position. ( if valid )

        :param seq_file:
        :return snp_seq:
        """

        # Read the sequence and store the data
        sequence = ""
        with open(seq_file) as file:
            for line in file:
                sequence += line.strip()

        # SNP len can only be 1
        if len(self.snp) > 1 or len(self.snp) < 0:
            raise Exception("\nSNP should be one char long")
        # The SNP must be a nucleotide
        if self.snp not in "ACTG":
            raise Exception("\nSNP should be a nucleotide: A,C,T or G")
        if self.pos in range(0, len(sequence)):

            # Put the snp in the correct position
            snp_seq = sequence[:self.pos] + self.snp + sequence[self.pos + 1:]
        else:
            raise Exception(F"\nThis index does not exist because the sequence is length: "
                            F"{(len(sequence) - 1)}")

        return snp_seq

    def dna_to_protein(self, sequence):

        """
        Function which translates DNA to Protein.

        :param sequence:
        :return protein_seq:
        """
        protein_seq = ""

        # The all-known codon table
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
            'TAC': 'Y', 'TAT': 'Y', 'TAA': '-', 'TAG': '-',
            'TGC': 'C', 'TGT': 'C', 'TGA': '-', 'TGG': 'W',
        }

        for i in range(0, len(sequence), 3):
            codon = sequence[i:i + 3]
            protein_seq += table[codon]
        print("Using protein sequence : \n{}\n".format(protein_seq))
        return protein_seq

    def get_align_score(self, alignment, snp_seq, column_len):
        """

        Function that takes an alignment and a sequence with an implemented SNP,
        and calculates the score based on preservation.

        :param alignment:
        :param snp_seq:
        :param column_len:
        :return:
        """
        scores = []
        aligner = Align.PairwiseAligner()
        for element in range(0, len(snp_seq)):
            alignments = aligner.align(alignment[:, element],
                                       snp_seq[element] * column_len)
            scores.append(alignments.score)
        return scores

    def write_results(self, score, column_len):
        """

        Function that takes scores and a position
        and writes the corresponding message.

        :param score:
        :param column_len:
        :return:
        """

        # Because the given position was specific for DNA, the pos needs to be corrected:
        # Example: position 3 of the DNA sequence is in protein 1, so 2 needs to be floored to 1.
        corrected_pos = math.floor(self.pos / 3)

        message = "Your score is not valid!"
        if score[corrected_pos] <= 0.1 * column_len:

            # Lower or equal to 10% similarity
            message = "So this is a pretty bad SNP, it has a lot of consequences!\n" \
                      "The AA might be pretty preserved!"
        if 0.1 * column_len < score[corrected_pos] <= 0.4 * column_len:

            # Between 10- and 40% similarity
            message = "So this SNP has some bad affects, but is not terrible."
        if 0.4 * column_len < score[corrected_pos] <= 0.8 * column_len:

            # Between 40- and 80% similarity
            message = "So the SNP has very few bad effect."
        if score[corrected_pos] > 0.8 * column_len:

            # Higher than 80% similarity
            message = "So this SNP has a neutral effect!"
        print(F"Severity of SNP at pos: {self.pos}, has score:"
              F" {score[corrected_pos]} / {column_len}.\n{message}")


def main(args):
    """
    The main function running all other functions
    :param args:
    :return:
    """
    print("\n- Starting Program -\n")

    # create a SNP Annotation object
    snp_obj = SnpAnnotation(args.in_file, args.SNP_Pos, args.SNP)

    # Pipeline:
    alignment, column_len = snp_obj.read_msa()  # Read MSA
    snp_seq = snp_obj.process_sequence(args.Sequence)  # Process the sequence to compare with
    protein_seq = snp_obj.dna_to_protein(snp_seq)  # From DNA --> Protein
    score = snp_obj.get_align_score(alignment, protein_seq, column_len)  # Get the alignment score
    snp_obj.write_results(score, column_len)  # Interpret results, and write them

    print("\n- End of Program -")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process a MSA.')
    parser.add_argument('-i', type=str, dest='in_file',
                        help='Please hand over a .msf file')
    parser.add_argument('-seq', type=str, dest='Sequence',
                        help='Please hand over a DNA or protein sequence file.')
    parser.add_argument('-pos', type=int, dest='SNP_Pos',
                        help='Please hand over a SNP position in the sequence.')
    parser.add_argument('-snp', type=str, dest='SNP',
                        help='Please hand over a valid SNP.')
    arguments = parser.parse_args()
    sys.exit(main(arguments))
