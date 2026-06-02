# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 08:52:48 2023

@author: yfu
"""

import sys, os
from utils import extract_prot_info_from_FASTAs

print("\n  Parsing FASTAs...")

# make directory for outputs
saveDir = os.path.join(os.getcwd(), sys.argv[1])
os.makedirs(saveDir, exist_ok=True)

# get fasta files
fasta_file = []
n = len(sys.argv)
for i in range(2, n):
    fasta_file.append(sys.argv[i])

# extract protein annotation
prot_info = extract_prot_info_from_FASTAs(fasta_file, out_file = os.path.join(saveDir, sys.argv[1] + '.csv'))

print("\n  Finished!")