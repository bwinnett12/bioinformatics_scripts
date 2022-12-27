from Bio import SeqIO
import sys
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def dna_to_rna(seq):
    dna_to_rna_table = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'U', '\n': ''}

    seq_to_return = ""
    for i in range(len(seq)):
        seq_to_return += dna_to_rna_table[seq[i]]

    return seq_to_return

def convert_dna_to_rna(file, out_file=""):
    my_records = []
    with open(out_file if out_file else "{path}_rna.fasta".format(path=os.path.splitext(file)[0]), "w") as f:
        for record in SeqIO.parse(file, "fasta"):
            record.description = ' '.join(record.description.split()[1:])
            trans_rec = SeqRecord(Seq(dna_to_rna(record.seq)), id=record.id, description=record.description)
            my_records.append(trans_rec)

        SeqIO.write(my_records, f, "fasta")

if __name__ == '__main__':
    if len(sys.argv) > 2:
        convert_dna_to_rna(sys.argv[1], out_file=sys.argv[2])
    else:
        convert_dna_to_rna(sys.argv[1])