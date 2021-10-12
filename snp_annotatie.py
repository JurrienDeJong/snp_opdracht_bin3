#!/bin/usr/python 3

"""
69
Nice
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

    def __init__(self, file, seq, pos, snp):
        self.file = file
        self.seq = seq
        self.pos = pos
        self.snp = snp

    def read_file(self):
        """

        Function which reads a .msf file with

        :param file_name:
        :return:
        """
        if self.file.endswith(".msf"):
            alignment = AlignIO.read(open(self.file), "msf")
            return alignment
        else:
            raise Exception("An incorrect filetype is given."
                            " Please hand over a .msf file")

    def dna_to_protein(self):
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
            'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
            'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
        }

        # Check if it is a DNA seq
        if not re.match("[BDEFHIJKLMNOPQRSUVWXYZ]", self.seq):
            for i in range(0, len(self.seq), 3):
                codon = self.seq[i:i + 3]
                protein_seq += table[codon]
            print("Using protein sequence : {}".format(protein_seq))
            return protein_seq
        else:
            # If it is a protein seq, return it.
            return self.seq

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

    def get_align_score(self, alignment, snp_seq):
        """

        Function that takes an alignment and a sequence with an implemented,
        and calculates the score based on preservation.

        :param alignment:
        :param snp_seq:
        :return:
        """
        scores = []
        aligner = Align.PairwiseAligner()
        for x in range(0, len(snp_seq)):
            alignments = aligner.align(alignment[:, x], snp_seq)
            scores.append(round(((alignments.score / len(snp_seq)) * 100), 2))
        return scores

    def write_results(self, percentage):
        """

        Function that takes percentage scores and a position
        and writes the corresponding message.

        :param percentage:
        :param pos:
        :return:
        """
        message = "Your score is not valid!"
        if percentage[self.pos] <= 10:
            message = "So this is a pretty bad SNP, it has a lot of consequences!\n" \
                      "The AA might be pretty preserved!"
        if 10 < percentage[self.pos] <= 40:
            message = "So this SNP has some bad affects, but is not terrible."
        if 40 < percentage[self.pos] <= 70:
            message = "So the SNP has very few bad effect."
        if percentage[self.pos] > 70:
            message = "So this SNP has a neutral effect! :)"
        print("Severity of SNP at pos: {}, has percentage score: {} %.\n{}"
              .format(self.pos, percentage[self.pos], message))


def main(arguments):
    """
    The main function running all other functions
    :param arguments:
    :return:
    """
    print("\n\n\n- Starting Program -\n")
    x = SnpAnnotation(arguments.in_file, arguments.Sequence, arguments.SNP_Pos, arguments.SNP)
    alignment = x.read_file()
    x.dna_to_protein()
    snp_seq = x.implement_snp()
    score = x.get_align_score(alignment, snp_seq)
    x.write_results(score)
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
