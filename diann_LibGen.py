#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 14:15:57 2023

@author: yfu
"""

import sys, os, shutil
from utils import getParams

# applying for a computing node
#os.system('bsub -P DIANN -n 16 -q standard -R "rusage[mem=8G]" -Is bash')

# read parameter file
paramFile = sys.argv[1]
params = getParams(paramFile)

# make a directory for output
saveDir = os.path.join(os.getcwd(), params["library_name"])
os.makedirs(saveDir, exist_ok=True)
# copy parameter file
shutil.copy(paramFile, os.path.join(saveDir,'diann_LibGen.params'))

# in silico digestion
digestion_rule = {'Trypsin_P': 'K*,R*',
                  'Trypsin': 'K*,R*,!*P',
                  'Lys_C': 'K*',
                  'Chymotrypsin': 'F*,Y*,W*,M*,L*,!*P',
                  'AspN': '*D',
                  'GluC': 'E*',}

# input fastas
fastas = params['fasta_file'].split(";")
fasta_cmd = ""
for fasta in fastas:
    temp = " --fasta " + fasta
    fasta_cmd = fasta_cmd + temp

# output file name
outputLib = os.path.join(saveDir, params['library_name']+".tsv")

# generate command
cmd_prt1 = (" --threads " + params['n_CPU_threads'] + " --verbose " + params['log_level'] +
           " --gen-spec-lib --predictor --fasta-search" + fasta_cmd + " --out-lib " + outputLib +
           " --cut " + digestion_rule[params['protease']] + " --missed-cleavages " + params['n_missed_cleavages'] +
           " --min-pep-len " + params['min_pep_length'] + " --max-pep-len " + params['max_pep_length'] +
           " --min-pr-mz " + params['min_prec_mz'] + " --max-pr-mz " + params['max_prec_mz'] + 
           " --min-pr-charge " + params['min_prec_charge'] + " --max-pr-charge " + params['max_prec_charge'] +
           " --min-fr-mz " + params['min_frag_mz'] + " --max-fr-mz " + params['max_frag_mz'])


cmd_prt2 = "" # for modifications
# fixed modifications
if params['C_carbamidomethylation'] == '1':
    cmd_prt2 = cmd_prt2 + " --unimod4"
# variable modifications
if params['M_oxidation'] == '1':
    cmd_prt2 = cmd_prt2 + " --var-mod UniMod:35,15.994915,M"
if params['N_term_M_excision'] == '1':
    cmd_prt2 = cmd_prt2 + " --met-excision"
# PTM analysis
if params['N_term_acetylation'] == '1':
    cmd_prt2 = cmd_prt2 + " --var-mod UniMod:1,42.010565,*n"
if params['Phosphorylation'] == '1':
    cmd_prt2 = cmd_prt2 + " --var-mod UniMod:21,79.966331,STY"
if params['K_GG_addcut'] == '1':
    cmd_prt2 = cmd_prt2 + " --var-mod UniMod:121,114.042927,K --no-cut-after-mod UniMod:121"

cmd_prt2 = cmd_prt2 + " --var-mods " + params['max_n_variable_mod']

cmd_prt3 = "" # for others

if params['additional_options'].lower() != 'null' and params['additional_options'] != "":# null (default); others, e.g. --var-mod Dimethyl,28.0313,KR
    cmd_prt3 = cmd_prt3 + " " + params['additional_options']

cmd = cmd_prt1 + cmd_prt2 + cmd_prt3

print('\n')
print(f'RUN:\n'+'apptainer exec /hpcf/authorized_apps/proteomics_apps/DIANN_CMD/DIA-NN_2_Linux_v1.9.2_v2.0.2/diann2.2.0_with_dotnet.sif /diann-2.2.0/diann-linux'+cmd)
os.system('apptainer exec /hpcf/authorized_apps/proteomics_apps/DIANN_CMD/DIA-NN_2_Linux_v1.9.2_v2.0.2/diann2.2.0_with_dotnet.sif /diann-2.2.0/diann-linux'+cmd)

