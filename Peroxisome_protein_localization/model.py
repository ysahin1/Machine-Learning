from Bio import motifs
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.svm import LinearSVC
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import ShuffleSplit
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

records = pd.read_csv("Reumann_Supplemental_Dataset_1A_Final.csv")
print(records["C-terminal_14_aa"].head())

classes = records["class"]

print ("hello world"1