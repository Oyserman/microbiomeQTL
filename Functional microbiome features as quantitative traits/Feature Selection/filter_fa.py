#!/usr/bin/env python

import sys
#from sets import Set
from Bio import SeqIO

def list_ids():
    """
    Return a set containing the identifiers presented in a file,
    line by line
    """

    # read the first file given and generate a set (faster iteration respect lists

    identifiers = set()

    with open(sys.argv[1], 'r') as fi:
        for line in fi:
            line = line.strip()
            identifiers.add(str(line))

    return identifiers

def filter():
    """
    Writes a file containing only the sequences with identifiers NOT
    present in a set of identifiers
    """

    identifiers = list_ids()

    with open(sys.argv[2]) as original_fasta, open(sys.argv[3], 'w') as corrected_fasta:
        records = SeqIO.parse(original_fasta, 'fasta')
        for record in records:
            print(record.id)
            if record.id not in identifiers:
                SeqIO.write(record, corrected_fasta, 'fasta')

if __name__ == '__main__':
    filter()
