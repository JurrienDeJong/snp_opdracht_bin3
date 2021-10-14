#!/bin/usr/python 3

"""
t.b.a.
"""

__author__ = "Jurrien de Jong"
__date__ = "8/10/2021"
__version__ = "V1.0"

import re
import sys
import argparse
from Bio import AlignIO
from Bio import Align


class SnpAnnotation:

    def __init__(self, file, pos, snp):
        self.file = file
        self.pos = pos
        self.snp = snp

    def read_data(self):
        """

        Function which reads a .msf file

        :param file_name:
        :return:
        """
        if self.file.endswith(".msf"):
            alignment = AlignIO.read(open(self.file), "msf")

            # Get the length of the column
            column_len = len(alignment[:, 1])
            return alignment, column_len
        else:
            raise Exception("An incorrect filetype is given."
                            " Please hand over a .msf file")

    def dna_to_protein(self, sequence):
        """
        Function which translates DNA to Protein if a DNA
        sequence is given.

        :param sequence:
        :return:
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
        print("Using protein sequence : {}".format(protein_seq))
        return protein_seq

    def implement_snp(self):
        """
        Function which implements a SNP in a given seq,
        at a given position. ( if valid )

        :param seq:
        :param pos:
        :param snp:
        :return:
        """
        if len(self.snp) > 1:
            raise Exception("SNP should be one char long")
        if self.pos in range(0, len(self.seq)):
            sequence = self.seq[:self.pos] + self.snp + self.seq[self.pos + 1:]
        else:
            raise Exception("This index does not exist"
                            " because the sequence is "
                            "length: {}".format(len(self.seq) - 1))
        return sequence

    def get_align_score(self, alignment, snp_seq, column_len):
        """

        Function that takes an alignment and a sequence with an implemented,
        and calculates the score based on preservation.

        :param alignment:
        :param snp_seq:
        :param column_len:
        :return:
        """
        scores = []
        aligner = Align.PairwiseAligner()
        for x in range(0, len(snp_seq)):
            alignments = aligner.align(alignment[:, x], snp_seq[x] * column_len)
            scores.append(alignments.score)
        return scores

    def write_results(self, score, column_len):
        """

        Function that takes percentage scores and a position
        and writes the corresponding message.

        :param score:
        :param column_len:
        :return:
        """
        message = "Your score is not valid!"
        if score[self.pos] <= 0.1 * column_len:

            # Lower or equal to 10% similarity
            message = "So this is a pretty bad SNP, it has a lot of consequences!\n" \
                      "The AA might be pretty preserved!"
        if 0.1 * column_len < score[self.pos] <= 0.4 * column_len:

            # Between 10- and 40% similarity
            message = "So this SNP has some bad affects, but is not terrible."
        if 0.4 * column_len < score[self.pos] <= 0.8 * column_len:

            # Between 40- and 80% similarity
            message = "So the SNP has very few bad effect."
        if score[self.pos] > 0.8 * column_len:

            # Higher than 80% similarity
            message = "So this SNP has a neutral effect! :)"
        print("Severity of SNP at pos: {}, has score: {} / {}.\n{}"
              .format(self.pos, score[self.pos], column_len, message))


def main(arguments):
    """
    The main function running all other functions
    :param arguments:
    :return:
    """
    print("\n- Starting Program -\n")

    # create a SNP Annotation object
    x = SnpAnnotation(arguments.in_file, arguments.SNP_Pos, arguments.SNP)

    alignment, column_len = x.read_data()
    x.dna_to_protein(sequence)
    snp_seq = x.implement_snp()
    score = x.get_align_score(alignment, snp_seq, column_len)
    x.write_results(score, column_len)
    print("\n- End of Program -")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process a MSA.')
    parser.add_argument('-i', type=str, dest='in_file',
                        help='Please hand over a .msf file')
    parser.add_argument('-seq', type=str, dest='Sequence',
                        help='Please hand over a DNA or protein sequence.')
    parser.add_argument('-pos', type=int, dest='SNP_Pos',
                        help='Please hand over a SNP position in the sequence.')
    parser.add_argument('-snp', type=str, dest='SNP',
                        help='Please hand over a valid SNP.')
    arguments = parser.parse_args()
    sys.exit(main(arguments))
