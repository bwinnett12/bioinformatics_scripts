from Bio import SeqIO
import sys


def find_gc(file):
    """
    Finds the GC content of each strand

    :param file: file to read gc_content
    :return: dict of the gc content of each entry
    """
    gc_dict = {}
    for record in SeqIO.parse(file, "fasta"):
        gc = 0
        for i in record.seq:
            if i == "C" or i == 'G':
                gc += 1

        gc_dict[record.id] = gc / len(record.seq)

    return gc_dict


def print_content(file):
    """
    Prints gc content of each entry

    :param file: file to print gc_content
    """
    for record in find_gc(file).items():
        print("{seq_id}: {seq_gc}".format(seq_id=record[0], seq_gc=record[1]))


if __name__ == '__main__':
    print_content(sys.argv[1])
