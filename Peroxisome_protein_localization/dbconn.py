import sqlite3
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
conn = sqlite3.connect('pssmdb.sqlite')

cur = conn.cursor()

cur.executescript('''
CREATE TABLE IF NOT EXISTS Pssm (
    id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
    protein_name TEXT UNIQUE,
    last_15_AA TEXT UNIQUE,
    pssm_table TEXT,
    sign INTEGER
);

CREATE TABLE IF NOT EXISTS Pssm (
    id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
    sign_no INTEGER UNIQUE, 
    title TEXT UNIQUE
);
''')

while True:
    fname = input('please type a file name: ')
    if not fname.endswith('.fasta'): 
        print ('please upload a fasta file')
        break
    try:
        reads = SeqIO.parse(fname, 'fasta')
    except:
        print ('The file is not a fasta')
        break
        
    def clean(line):
    # remove extra spaces
        t_string = line.rstrip()
        t_string = re.sub(' +',' ',t_string)
        t_string = t_string.split(' ')
        return t_string

    def normalize(pssm_matrix):
    # normalize pssm matrix
    # it uses 1/1+e**(-x) normalization formula
        norm_rows = []
        for i in pssm_matrix:
            norm_row = []
            for ii in i:
                norm = float(int(1)/(int(1)+(math.exp(int(ii)*int(-1)))))
                norm_row.append(norm)
            norm_rows.append(norm_row)
        return norm_rows

    def write_ind_rec(record):
        # write sequences individually 
        output_file = "individiual.fasta"
        output_handle = open(output_file, "w")
        SeqIO.write(record, output_handle, "fasta")
        output_handle.close()
        return print ("Fasta file has been generated : {}".format(record.id))

    all_pssm_matrix = []

    def check_if_file_exist(record_id):
        pssm_file = "./pssm/pssm_{}.csv".format(record_id)
        if os.path.exists(pssm_file):
            aa_sequence = []
            pssm_matrix = []
            with open(pssm_file,'r') as pssm_file:
                t_file_lines = pssm_file.readlines()
                for line in t_file_lines:
                    t_string = clean(line)
                    if len(t_string) >= 44:
                        #print (t_string)
                        aa_sequence.append(t_string[2])
                        pssm_matrix.append(t_string[3:3+20])
            pssm = pssm_matrix[35:]
            pssm_matrix = normalize(pssm)
            all_pssm_matrix.append(pssm_matrix)
            os.remove("./pssm/pssm_{}.csv".format(record_id))
            print("File Removed!")
        else:
            no_hits.append(record_id)
            return print ("File does not exist: {} \n".format(record_id))
        
    def get_psi_results(reads):
        number_sequence = 0
        for record in reads:
            number_sequence += 1
            #write a sequence into seperate file
            write_ind_rec(record)
            # send the protein sequence to psi-blast
            psi_blast_cmd = "psiblast \
            -db /media/yunus/TOSHIBA1TB/nr/nr \
            -evalue 0.001 \
            -query /home/yunus/Desktop/MEGA/ml_p/pwm/individiual.fasta \
            -out_ascii_pssm pssm/pssm_{}.csv \
            -num_iterations 3"
            psi_blast_cmd = psi_blast_cmd.format(record.id)
            p = subprocess.Popen(psi_blast_cmd, stdout=subprocess.PIPE, shell=True)
            output, err = p.communicate()  
            #This makes the wait possible
            p_status = p.wait()
            #check if file exist
            check_if_file_exist(record.id)
        return (all_pssm_matrix, no_hits, number_sequence)

    all_pssm_matrix, no_blast, number_sequence = get_psi_results(reads)
    print ("\n")
    print ("The number of processed sequence: {}".format(number_sequence))
    print ("\n")
    print ("Length of the generated array: {}".format(len(all_pssm_matrix)))
    print ("\n")
    print ("The number of no psi-Blast: {}".format(len(no_blast)))