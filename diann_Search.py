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
saveDir = os.path.join(os.getcwd(), params["outputDir"])
os.makedirs(saveDir, exist_ok=True)
# copy parameter file
shutil.copy(paramFile, os.path.join(saveDir,'diann_Search.params'))

# input fastas
input_runs = params['input_runs']

iuput_runs_cmd = ""
for run in input_runs:
    temp = " --f " + run
    iuput_runs_cmd = iuput_runs_cmd + temp

fastas = params['fasta_file'].split(";")
fasta_cmd = ""
for fasta in fastas:
    temp = " --fasta " + fasta
    fasta_cmd = fasta_cmd + temp

# output file name
outputReport = os.path.join(saveDir, "report.parquet")
outputLib = os.path.join(saveDir, "report-lib.parquet")

# generate command
cmd_prt0 = (" --threads " + params['n_CPU_threads'] + " --verbose " + params['log_level'] + 
            iuput_runs_cmd + fasta_cmd + " --qvalue " + params['fdr_prec'] + " --lib " + params['SpecLib'] + 
           " --gen-spec-lib" + " --out-lib " + outputLib + " --out " + outputReport)

cmd_prt1 = "" # for modifications

cmd_prt1 = cmd_prt1 + " --var-mods " + params['max_n_variable_mod']

# fixed modifications
if params['C_carbamidomethylation'] == '1':
    cmd_prt1 = cmd_prt1 + " --unimod4"
# variable modifications
if params['M_oxidation'] == '1':
    cmd_prt1 = cmd_prt1 + " --var-mod UniMod:35,15.994915,M"
if params['N_term_M_excision'] == '1':
    cmd_prt1 = cmd_prt1 + " --met-excision"
# PTM analysis
if params['N_term_acetylation'] == '1':
    cmd_prt1 = cmd_prt1 + " --var-mod UniMod:1,42.010565,*n"
if params['Phosphorylation'] == '1':
    cmd_prt1 = cmd_prt1 + " --var-mod UniMod:21,79.966331,STY"
if params['K_GG_addcut'] == '1':
    cmd_prt1 = cmd_prt1 + " --var-mod UniMod:121,114.042927,K --no-cut-after-mod UniMod:121"

cmd_prt2 = "" # Identification related parameters 
# if use MBR
if params['MBR'] == '1':
    cmd_prt2 = cmd_prt2 + " --reanalyse"
if params['ms2_tolerance'] != '0': # if automatic inference is not chosen
    cmd_prt2 = cmd_prt2 + " --mass-acc " + params['ms2_tolerance']
if params['ms1_tolerance'] != '0': # if automatic inference is not chosen
    cmd_prt2 = cmd_prt2 + " --mass-acc-ms1 " + params['ms1_tolerance']
if params['n_scans_per_peak'] != '0': # if automatic inference is not chosen
    cmd_prt2 = cmd_prt2 + " --window " + params['n_scans_per_peak']
# whether to use precursor isotopes
# if params['use_isotopologues'] == '0':
    # cmd_prt2 = cmd_prt2 + " --no-isotopes"
# # whether to remove interference precursor ID
# if params['noise_precursor_ID_removal'] == '0':
    # cmd_prt2 = cmd_prt2 + " --int-removal 0"
# if unrelated runs
if params['unrelated_runs'] == '1':
    cmd_prt2 = cmd_prt2 + " --individual-mass-acc --individual-windows"
# if report lib info
# if params['report_LibInfo'] == '1':
    # cmd_prt2 = cmd_prt2 + " --report-lib-info"
if params['scoring'] == '2':# 1 = Generic (default); 2 = Peptidoforms; 3 = Proteoforms
    cmd_prt2 = cmd_prt2 + " --peptidoforms"
elif params['scoring'] == '3':
    cmd_prt2 = cmd_prt2 + " --proteoforms"
if params['machine_learning'] == '1':# 1 = Linear classifier; 2 = NNs (fast); 3 = NNs (cross-validated) (default)
    cmd_prt2 = cmd_prt2 + " --no-nn"
elif params['machine_learning'] == '2':
    cmd_prt2 = cmd_prt2 + " --fast-ml"

cmd_prt3 = "" # for quantification and protein inference

# protein inference

# if params['heuristic_prot_inference'] == '1':
    # cmd_prt3 = cmd_prt3 + " --relaxed-prot-inf"
    
if params['prot_inference'] == '0':
    cmd_prt3 = cmd_prt3 + " --no-prot-inf"
if params['proteotypicity'] == '1':# 1 = Isofrom IDs; 2 = Protein names (from FASTA); 3 = Genes (Species-specific); 4 = Genes (default)
    cmd_prt3 = cmd_prt3 + " --pg-level 0"
elif params['proteotypicity'] == '2':
    cmd_prt3 = cmd_prt3 + " --pg-level 1"
elif params['proteotypicity'] == '3':
    cmd_prt3 = cmd_prt3 + " --pg-level 2 --species-genes"

# quantification

if params['quantification_strategy'] == '1':# 1 = Legacy (direct); 2 = QuantUMS (high accuracy); 3 = QuantUMS (high precision) (default)
    cmd_prt3 = cmd_prt3 + " --direct-quant"
elif params['quantification_strategy'] == '2':
    cmd_prt3 = cmd_prt3 + " --high-acc"
# if params['use_peak_hight'] == '1':
    # cmd_prt3 = cmd_prt3 + " --peak-height"
# elif params['use_peak_hight'] == '0':
    # if params['robust_LC'] == '1':
        # cmd_prt3 = cmd_prt3 + " --peak-center"
    # if params['frag_noise_removal'] == '0':
        # cmd_prt3 = cmd_prt3 + " --no-ifs-removal"

# normalization

if params['cross_run_normalization'] == '1':# 1 = Global; 2 = RT-dependent (default); 3 = RT & signal-dep. (experimental); 4 = off
    cmd_prt3 = cmd_prt3 + " --global-norm"
elif params['cross_run_normalization'] == '3':
    cmd_prt3 = cmd_prt3 + " --sig-norm"
elif params['cross_run_normalization'] == '4':
    cmd_prt3 = cmd_prt3 + " --no-norm"
# if params['cross_run_normalization'] == '0':
    # cmd_prt3 = cmd_prt3 + " --no-norm"
# elif params['cross_run_normalization'] == '2':
    # cmd_prt3 = cmd_prt3 + " --global-norm"

cmd_prt4 = "" # for others

# use existing .quant files if available
if params['use_quant'] == '1':
    cmd_prt4 = cmd_prt4 + " --use-quant"

if params['re_annotate'] == '1':
    cmd_prt4 = cmd_prt4 + " --reannotate"

if params['quantities_matrices'] == '1':
    cmd_prt4 = cmd_prt4 + " --matrices"

if params['XICs'] == '1':
    cmd_prt4 = cmd_prt4 + " --xic"

if params['library_generation'] == '1':# 1 = IDs profiling; 2 = IDs, RT & IM profiling (default); 3 = Smart profiling; 4 = Full profiling
    cmd_prt4 = cmd_prt4 + " --id-profiling"
elif params['library_generation'] == '2':
    cmd_prt4 = cmd_prt4 + " --rt-profiling"
elif params['library_generation'] == '3':
    cmd_prt4 = cmd_prt4 + " --smart-profiling"

if params['speed_RAM_usage'] == '2':# 1 = Optimal results (default); 2 = Low RAM usage; 3 = Low RAM & high speed; 4 = Ultra-fast
    cmd_prt4 = cmd_prt4 + " --min-corr 1.0 --time-corr-only"
elif params['speed_RAM_usage'] == '3':
    cmd_prt4 = cmd_prt4 + " --min-corr 2.0 --time-corr-only"
elif params['speed_RAM_usage'] == '4':
    cmd_prt4 = cmd_prt4 + " --min-corr 2.0 --time-corr-only --extracted-ms1"

if params['additional_options'].lower() != 'null' and params['additional_options'] != "":# null (default); others, e.g. --var-mod Dimethyl,28.0313,KR
    cmd_prt4 = cmd_prt4 + " " + params['additional_options']

#cmd_prt4 = cmd_prt4 + " --xic"

cmd = cmd_prt0 + cmd_prt1 + cmd_prt2 + cmd_prt3 + cmd_prt4

print('\n')
print('cmd\napptainer exec /hpcf/authorized_apps/proteomics_apps/DIANN_CMD/DIA-NN_2_Linux_v1.9.2_v2.0.2/diann2.2.0_with_dotnet.sif /diann-2.2.0/diann-linux'+cmd)
os.system('apptainer exec /hpcf/authorized_apps/proteomics_apps/DIANN_CMD/DIA-NN_2_Linux_v1.9.2_v2.0.2/diann2.2.0_with_dotnet.sif /diann-2.2.0/diann-linux'+cmd)

