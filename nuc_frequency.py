from Bio import SeqIO
import sys


def countNucFrequency(seq):
    """
    Finds counts of each nucleotide. Extends past ACGT

    :param seq: sequence to count nucleotide frequencies
    :return: dict of nucleotide frequencies
    """
    tmpFreqDict = {"A": 0, "C": 0, "G": 0, "T": 0,
                   "U": 0, "W": 0, "S": 0, "M": 0,
                   "K": 0, "R": 0, "Y": 0, "B": 0,
                   "D": 0, "H": 0, "V": 0, "N": 0}

    for nuc in seq:
        tmpFreqDict[nuc] += 1
    return tmpFreqDict


def get_nuc_freq(file):
    """
    Finds the nucleotide frequencies of each strand in a fasta file

    :param file: a fasta file of nucleotide strands
    """
    for record in SeqIO.parse(file, "fasta"):
        print(record.description)
        print('\n'.join(["{base}: {numb}".format(base=key, numb=str(val))
                         for key, val in countNucFrequency(record.seq).items() if val > 0]))
        print(" ")


if __name__ == '__main__':
    get_nuc_freq(sys.argv[1])
