#generete pssm tables.
import math
import subprocess 
import re
import os
from Bio import SeqIO

class Generate_pssm:

    def __init__(self,fname):
        number_sequence = 0
        self.fname = fname

    reads = SeqIO.parse(self.fname, 'fasta')
    for record in reads:
        print (record.id)
deneme = Generate_pssm(fname='results.fasta')
