#!/biodata4/home/cliu/cpgview/bin/python

from Bio import GenBank
import sys

file_gb = sys.argv[1]
with open(file_gb) as handle:
    for record in GenBank.parse(handle):
        #print(record.accession)
        print(record.organism)
