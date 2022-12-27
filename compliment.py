from Bio import SeqIO
import sys
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def dna_compliment(seq):
    """
    Finds the compliment of a given sequence and returns it

    :param seq: The sequence to find the compliment of
    :return: String containing the sequence
    """
    dna_compliment_table = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '\n': ''}
    seq_to_return = ""

    for i in range(len(seq)):
        seq_to_return += dna_compliment_table[seq[i]]

    return seq_to_return[::-1]


def convert_compliment(file, out_file=""):
    """
    Finds the compliment of each strand of a fasta file

    :param file: input file to find compliments of
    :param out_file: [optional] name of output file
    """
    my_records = []
    with open(out_file if out_file else "{path}_compliment.fasta".format(path=os.path.splitext(file)[0]), "w") as f:
        for record in SeqIO.parse(file, "fasta"):
            record.description = ' '.join(record.description.split()[1:])
            record.description += "|compliment"
            trans_rec = SeqRecord(Seq(dna_compliment(record.seq)), id=record.id, description=record.description)
            my_records.append(trans_rec)

        SeqIO.write(my_records, f, "fasta")


if __name__ == '__main__':
    if len(sys.argv) > 2:
        convert_compliment(sys.argv[1], out_file=sys.argv[2])
    else:
        convert_compliment(sys.argv[1])
