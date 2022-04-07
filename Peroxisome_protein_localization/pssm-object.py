import math
import subprocess 
import re
import os
from Bio.Seq import Seq
import numpy as np
import pandas as pd
from numpy import asarray
from numpy import save
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

class Generate_pssm:
    all_pssm_matrix = list()
    no_hits = list()
    number_sequence = 0
    norm_rows = []
    def __init__(self,fname = "file.fasta"):
        self.fname = fname # assign a variable name to file
        reads = SeqIO.parse(self.fname, 'fasta') # parse fasta file
        for record in reads:
            self.number_sequence += 1
            self.record_id = record.id
            self.record = record
            self.parse(record)
            self.blast(self.record_id)
            self.check_if_file_exist(self.record_id)
            return None
    # split each sequence
    def blast(self, record_id):
      self.record_id = record_id
      psi_blast_cmd = "psiblast -db /media/yunus/TOSHIBA1TB/nr/nr -evalue 0.001 -query /home/yunus/Desktop/projects/Peroxisome_svm_model/record.fasta -out_ascii_pssm /home/yunus/Desktop/projects/Peroxisome_svm_model/pssm/pssm_{}.csv -num_iterations 3"  
      psi_blast_cmd = psi_blast_cmd.format(self.record_id)
      process = subprocess.Popen(psi_blast_cmd, stdout=subprocess.PIPE, shell=True)
      output, err = process.communicate()  
      #This makes the wait possible
      process_status = process.wait()
            #check if file exist
      return print ('csv file has been generated')

    def parse(self, record):
        OUTPUT_FILE = "record.fasta"
        output_handle = open(OUTPUT_FILE, "w") # YENİ OLUŞTURULAN DOSYAYI AÇ
        SeqIO.write(self.record, output_handle, "fasta") # TEK SEQUENCE'I YAZ
        output_handle.close()
        return print ("Fasta file has been generated : {}".format(record.id))
    def clean(self, line):
        # remove extra spaces
        self.line = line
        t_string = self.line.rstrip()
        t_string = re.sub(' +',' ',t_string)
        t_string = t_string.split(' ')
        return t_string
    def normalize(self, pssm_matrix):
        self.pssm_matrix = pssm_matrix
        # normalize pssm matrix
        # it uses 1/1+e**(-x) normalization formula
        for i in self.pssm_matrix:
            norm_row = []
            for ii in i:
                norm = float(int(1)/(int(1)+(math.exp(int(ii)*int(-1)))))
                norm_row.append(norm)
            self.norm_rows.append(norm_row)
        return self.norm_rows

    def check_if_file_exist(self, record_id):
        self.record_id = record_id
        pssm_file = "./pssm/pssm_{}.csv".format(self.record_id)
        if os.path.exists(pssm_file):
            self.aa_sequence = []
            self.pssm_matrix = []
            with open(pssm_file,'r') as pssm_file:
                t_file_lines = pssm_file.readlines()
                for line in t_file_lines:
                    self.t_string = self.clean(line)
                    if len(self.t_string) >= 44:
                        print (self.t_string)
                        self.aa_sequence.append(self.t_string[2])
                        self.pssm_matrix.append(self.t_string[3:3+20])
                        pssm = self.pssm_matrix[35:]
                        pssm_matrix = self.normalize(pssm)
                        self.all_pssm_matrix.append(pssm_matrix)
                        os.remove("./pssm/pssm_{}.csv".format(self.record_id))
                        print("File Removed!")
                        return self.all_pssm_matrix
                    else:
                       self.no_hits.append(self.record_id)
                       return print ("File does not exist: {} \n".format(self.record_id))
deneme = Generate_pssm('results.fasta')
deneme.all_pssm_matrix
