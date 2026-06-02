#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 17:06:25 2023

@author: yfu, zfyuan
"""

import sys, os, shutil
from utils import getParams, Tee, getLoadingBias, normalization, read_report, get_uniprot, get_KRT, reorder_dataframe
import numpy as np, pandas as pd

#R CMD ldd /home/yfu/.conda/envs/diann_iq/lib/R/library/stats/libs/stats.so
#os.environ['R_HOME'] = 'C:/Users\yfu\AppData\Local\miniconda3\envs\diann_iq\Lib/R'
#os.environ['R_USER'] = 'C:/Users\yfu\AppData\Local\miniconda3\envs\diann_iq\lib\site-packages/rpy2'

#import rpy2.situation
#for row in rpy2.situation.iter_info():
#     print(row)

# import rpy2.robjects as robj
# from rpy2.robjects.packages import importr

# read parameter file
paramFile = sys.argv[1]
# paramFile = "/research/rgs01/home/clusterHome/yfu/projects/PengLab/06052024_APOE_IP_human_brain/diann_iqQuan.params"
params = getParams(paramFile)

# make a directory for output
saveDir = os.path.join(os.getcwd(), params["outputDir"])
# saveDir = os.path.join(os.path.dirname(paramFile), "testQuan")
os.makedirs(saveDir, exist_ok=True)
# copy parameter file
shutil.copy(paramFile, os.path.join(saveDir,'diann_iqQuan.params'))

# make the intermediate directory 
interDir = saveDir+'/intermediate'
os.makedirs(interDir, exist_ok=True)
# make the publication directory 
pubDir = saveDir+'/publication'
os.makedirs(pubDir, exist_ok=True)

#############################
# iq (r package) processing #
#############################

with Tee(os.path.join(saveDir,'log.txt')):
    
    # read DIA-NN report.tsv file
    diann_report = params["diann_report"]
    
    print("\nReading DIA-NN report file: " + diann_report)
    # diann_report_df = pd.read_table(diann_report)
    diann_report_df = read_report(diann_report)
    if "Run" not in list(diann_report_df.columns) or len(diann_report_df["Run"])==0:
        sys.exit( "\nno report.parquet/report.tsv in the path: {}".format(os.path.split(diann_report)[0]) )
    params_Qvalue = {key: value for key, value in params.items() if 'Q.Value' in key}
    
    # '''
    if params['QuanMethod'] == '1':
    # '''
    # if params['QuanMethod'] == '1':
        print("\nUsing PG.MaxLFQ values from DIA-NN report...")
        for key, value in params_Qvalue.items():
            diann_report_df = diann_report_df[diann_report_df[key] < float(value)]
        extracted_cols = ['Run', 'Protein.Group', #"Protein.Ids", 
                          'Protein.Names', 'Genes', 'PG.MaxLFQ']
        all_prot_quan_long = diann_report_df.loc[:, extracted_cols].drop_duplicates(ignore_index=True)
        all_prot_quan = all_prot_quan_long.pivot(index=['Protein.Group', 'Protein.Names','Genes'], 
                                                 columns='Run', values='PG.MaxLFQ')
        all_prot_quan.reset_index(drop=False, inplace=True)
        all_prot_quan.iloc[:, 3:] = np.log2(all_prot_quan.iloc[:, 3:] + 1)
    elif params['QuanMethod'] == '0':
        print("\nUsing PG.TopN values from DIA-NN report...")
        for key, value in params_Qvalue.items():
            diann_report_df = diann_report_df[diann_report_df[key] < float(value)]
        extracted_cols = ['Run', 'Protein.Group', #"Protein.Ids", 
                          'Protein.Names', 'Genes', 'PG.TopN']
        all_prot_quan_long = diann_report_df.loc[:, extracted_cols].drop_duplicates(ignore_index=True)
        all_prot_quan = all_prot_quan_long.pivot(index=['Protein.Group', 'Protein.Names','Genes'], 
                                                 columns='Run', values='PG.TopN')
        all_prot_quan.reset_index(drop=False, inplace=True)
        all_prot_quan.iloc[:, 3:] = np.log2(all_prot_quan.iloc[:, 3:] + 1)
    else:
        sys.exit("\nQuantification method not supported, please check your parameters!")
    
    # ----precur----
    extracted_cols = ['Run', "Precursor.Id", 'Protein.Group', 
                      'Protein.Names', 'Genes', 'Precursor.Quantity']
    all_prec_quan_long = diann_report_df.loc[:, extracted_cols].drop_duplicates(ignore_index=True)
    all_prec_quan = all_prec_quan_long.pivot(index=['Precursor.Id', 'Protein.Group', 'Protein.Names','Genes'], 
                                             columns='Run', values='Precursor.Quantity')
    all_prec_quan.reset_index(drop=False, inplace=True)
    all_prec_quan.iloc[:, 4:] = np.log2(all_prec_quan.iloc[:, 4:] + 1)
    
    extracted_cols = ['Run', "Precursor.Id", 'Protein.Group', 
                      'Protein.Names', 'Genes', 'RT']
    all_prec_rt_long = diann_report_df.loc[:, extracted_cols].drop_duplicates(ignore_index=True)
    all_prec_rt = all_prec_rt_long.pivot(index=['Precursor.Id', 'Protein.Group', 'Protein.Names','Genes'], 
                                             columns='Run', values='RT')
    all_prec_rt.reset_index(drop=False, inplace=True)
    # ----precur----
    
    # get sample names
    sample_names = all_prot_quan.columns.values.tolist()[3:]
    # sample_names_new = list(map(lambda x: os.path.basename(x).split(".")[0], sample_names))
    sample_names_new = list(map(lambda x: os.path.basename(x).split(".")[0].replace("-", "_"), sample_names))
    all_prot_quan.rename(columns = dict(zip(sample_names, sample_names_new)), inplace = True)
    all_prot_quan[sample_names_new] = all_prot_quan[sample_names_new].replace(0, np.nan)# replace 0 with np.nan for sample columns
    n_samples = len(sample_names_new)
    # ----precur----
    all_prec_quan.rename(columns = dict(zip(sample_names, sample_names_new)), inplace = True)
    all_prec_quan[sample_names_new] = all_prec_quan[sample_names_new].replace(0, np.nan)# replace 0 with np.nan for sample columns
    
    all_prec_rt.rename(columns = dict(zip(sample_names, sample_names_new)), inplace = True)
    all_prec_rt[sample_names_new] = all_prec_rt[sample_names_new].replace(0, np.nan)# replace 0 with np.nan for sample columns
    # ----precur----

    # read protein annotation file
    if params["prot_anno_file"] == "0":
        print("\nNo protein annotation file found...")
        all_prot_quan_anno = all_prot_quan
        # ----precur----
        all_prec_quan_anno = all_prec_quan
        
        all_prec_rt_anno = all_prec_rt
        # ----precur----
    else:
        print("\nAdd protein annotation...")
        print("  using protein annotation file: " + params["prot_anno_file"])
        prot_info = pd.read_csv(params["prot_anno_file"])
        prot_info.drop(columns=['Protein.Names','Database'], inplace = True)
        prot_info.rename(columns={"Protein.Ids": "Protein.Group.Rep"}, inplace = True)
        prot_info = prot_info.drop_duplicates(subset=['Protein.Group.Rep'], ignore_index=True)
        # add protein annotation to quantification matrix
        all_prot_quan.insert(1, 'Protein.Group.Rep', all_prot_quan['Protein.Group'].str.split(pat=';', n=1, expand=True)[0])
        all_prot_quan_anno = prot_info.merge(all_prot_quan, how = 'right')
        all_prot_quan_anno.drop(columns=['Protein.Group.Rep'], inplace = True)
        # ----precur----
        all_prec_quan.insert(1, 'Protein.Group.Rep', all_prec_quan['Protein.Group'].str.split(pat=';', n=1, expand=True)[0])
        all_prec_quan_anno = prot_info.merge(all_prec_quan, how = 'right')
        all_prec_quan_anno.drop(columns=['Protein.Group.Rep'], inplace = True)
        all_prec_quan_anno = reorder_dataframe(all_prec_quan_anno,'Precursor.Id','Accession.Number',1)
        
        all_prec_rt.insert(1, 'Protein.Group.Rep', all_prec_rt['Protein.Group'].str.split(pat=';', n=1, expand=True)[0])
        all_prec_rt_anno = prot_info.merge(all_prec_rt, how = 'right')
        all_prec_rt_anno.drop(columns=['Protein.Group.Rep'], inplace = True)
        all_prec_rt_anno = reorder_dataframe(all_prec_rt_anno,'Precursor.Id','Accession.Number',1)
        # ----precur----

    
    print("\nExtract peptide information...")
    # peptide-level information
    extracted_cols = ['Protein.Group', 'Protein.Names', 'Genes', 'Stripped.Sequence']
    prot_pep_info = diann_report_df.loc[:, extracted_cols].drop_duplicates(ignore_index=True)
    if all_prot_quan_anno['Genes'].isnull().sum() != all_prot_quan_anno.shape[0]:
        all_prot_pep = prot_pep_info.groupby(['Protein.Group', 'Protein.Names', 'Genes'], as_index=False, dropna=False).agg({'Stripped.Sequence': ', '.join}).reset_index(drop=True)
        all_prot_pep.insert(4, 'n.Identified.Peptides', all_prot_pep['Stripped.Sequence'].str.count(",") + 1)
    elif all_prot_quan_anno['Genes'].isnull().sum() == all_prot_quan_anno.shape[0]:
        all_prot_pep = prot_pep_info.groupby(['Protein.Group', 'Protein.Names'], as_index=False, dropna=False).agg({'Stripped.Sequence': ', '.join}).reset_index(drop=True)
        all_prot_pep.insert(3, 'n.Identified.Peptides', all_prot_pep['Stripped.Sequence'].str.count(",") + 1)
    # add peptide information
    all_prot_quan_anno = all_prot_pep.merge(all_prot_quan_anno, how = 'right')
    all_prot_quan_anno.rename(columns={"Stripped.Sequence": "Identified.Peptides"}, inplace = True)
    
    # Precursor.Mz information
    extracted_cols = ['Protein.Group', 'Protein.Names', 'Genes', 'Precursor.Mz']
    prot_pep_info = diann_report_df.loc[:, extracted_cols].drop_duplicates(ignore_index=True)
    if all_prot_quan_anno['Genes'].isnull().sum() != all_prot_quan_anno.shape[0]:
        all_prot_pep = prot_pep_info.groupby(['Protein.Group', 'Protein.Names', 'Genes'], as_index=False, dropna=False).agg({'Precursor.Mz': lambda x: ', '.join(map(str, x))}).reset_index(drop=True)
    elif all_prot_quan_anno['Genes'].isnull().sum() == all_prot_quan_anno.shape[0]:
        all_prot_pep = prot_pep_info.groupby(['Protein.Group', 'Protein.Names'], as_index=False, dropna=False).agg({'Precursor.Mz': lambda x: ', '.join(map(str, x))}).reset_index(drop=True)
    # add Precursor.Mz information
    all_prot_quan_anno = all_prot_pep.merge(all_prot_quan_anno, how = 'right')
    
    # ----reorder----
    all_prot_quan_anno = reorder_dataframe(all_prot_quan_anno,'Accession.Number','Genes')
    all_prec_quan_anno = reorder_dataframe(all_prec_quan_anno,'Protein.Names','Accession.Number')
    all_prec_quan_anno = reorder_dataframe(all_prec_quan_anno,'Genes','Protein.Description')
    
    all_prec_rt_anno = reorder_dataframe(all_prec_rt_anno,'Protein.Names','Accession.Number')
    all_prec_rt_anno = reorder_dataframe(all_prec_rt_anno,'Genes','Protein.Description')
    # ----reorder----
    
    all_prot_quan_anno_csv = pd.concat([all_prot_quan_anno.drop(columns=sample_names_new), all_prot_quan_anno.loc[:, sample_names_new].apply(np.exp2) - 1], axis=1)
    all_prot_quan_anno_csv.to_csv(os.path.join(interDir, "all_prot_quan.csv"), index=False)
    # ----precur----
    all_prec_quan_anno_csv = pd.concat([all_prec_quan_anno.drop(columns=sample_names_new), all_prec_quan_anno.loc[:, sample_names_new].apply(np.exp2) - 1], axis=1)
    all_prec_quan_anno_csv.to_csv(os.path.join(interDir, "all_prec_quan.csv"), index=False)
    
    all_prec_rt_anno_csv = all_prec_rt_anno# pd.concat([all_prec_rt_anno.drop(columns=sample_names_new), all_prec_rt_anno.loc[:, sample_names_new].apply(np.exp2) - 1], axis=1)
    all_prec_rt_anno_csv.to_csv(os.path.join(interDir, "all_prec_rt.csv"), index=False)
    # ----precur----
    
    # + extra steps +
    # + extra step1 +: Remove con|, cu|, Krt proteins;  Keep certain customized proteins (e.g. Abeta peptide and PSEN1)
    protein_include_list = ['cu|P05067-107|A4_Abeta40_HUMAN','cu|P05067-108|A4_Abeta42_HUMAN','cu|P05067-109|A4_Abeta40_','cu|P05067-110|A4_Abeta42_','cu|P49768-101|PSN1_5xFAD_HUMAN','cu|P49768-102|PSN1_deltaE9_HUMAN']
    #all_prot_quan_anno["ProtType"] = all_prot_quan_anno.apply(lambda x: get_uniprot(x["Accession.Number"], protein_include_list), axis=1)
    all_prot_quan_anno["ProtType"] = all_prot_quan_anno.apply(lambda x: get_uniprot(str(x["Accession.Number"]), protein_include_list) if pd.notna(x["Accession.Number"]) else "Unknown", axis=1)
    all_prot_quan_anno["KRT"] = all_prot_quan_anno.apply(lambda x: get_KRT(str(x["Genes"])), axis=1)
    all_prot_quan_anno = all_prot_quan_anno.copy(deep=False).loc[all_prot_quan_anno["ProtType"] == "uniprot"]
    all_prot_quan_anno = all_prot_quan_anno.copy(deep=False).loc[all_prot_quan_anno["KRT"] == "nonKRT"]
    all_prot_quan_anno.drop(columns=['ProtType','KRT'], inplace = True)
    # + extra steps +
    # ----precur----
    all_prec_quan_anno["ProtType"] = all_prec_quan_anno.apply(lambda x: get_uniprot(str(x["Accession.Number"]), protein_include_list) if pd.notna(x["Accession.Number"]) else "Unknown", axis=1)
    all_prec_quan_anno["KRT"] = all_prec_quan_anno.apply(lambda x: get_KRT(str(x["Genes"])), axis=1)
    all_prec_quan_anno = all_prec_quan_anno.copy(deep=False).loc[all_prec_quan_anno["ProtType"] == "uniprot"]
    all_prec_quan_anno = all_prec_quan_anno.copy(deep=False).loc[all_prec_quan_anno["KRT"] == "nonKRT"]
    all_prec_quan_anno.drop(columns=['ProtType','KRT'], inplace = True)
    
    all_prec_rt_anno["ProtType"] = all_prec_rt_anno.apply(lambda x: get_uniprot(str(x["Accession.Number"]), protein_include_list) if pd.notna(x["Accession.Number"]) else "Unknown", axis=1)
    all_prec_rt_anno["KRT"] = all_prec_rt_anno.apply(lambda x: get_KRT(str(x["Genes"])), axis=1)
    all_prec_rt_anno = all_prec_rt_anno.copy(deep=False).loc[all_prec_rt_anno["ProtType"] == "uniprot"]
    all_prec_rt_anno = all_prec_rt_anno.copy(deep=False).loc[all_prec_rt_anno["KRT"] == "nonKRT"]
    all_prec_rt_anno.drop(columns=['ProtType','KRT'], inplace = True)
	# ----precur----
    
    # # params["min_nonmiss_pct"], params["min_nonmiss_num"]
    # if "min_nonmiss_pct" not in params.keys():
        # params["min_nonmiss_pct"] = "0"
    # if "min_nonmiss_num" not in params.keys():
        # params["min_nonmiss_num"] = "3"
    # if n_samples-float(params["min_nonmiss_num"])<0:
        # params["min_nonmiss_num"] = "0"
    
    # remove missing values
    rm_index = all_prot_quan_anno.loc[:, sample_names_new].isnull().values.sum(axis=1) > min( n_samples * (float(params["allowed_miss_pct"]) / 100), n_samples-1 )
    # rm_index = all_prot_quan_anno.loc[:, sample_names_new].isnull().values.sum(axis=1) > min( n_samples * (1-float(params["min_nonmiss_pct"]) / 100), n_samples-float(params["min_nonmiss_num"]) )
    all_prot_quan_anno_rmNA = all_prot_quan_anno[~rm_index]
    # ----precur----
    rm_index = all_prec_quan_anno.loc[:, sample_names_new].isnull().values.sum(axis=1) > min( n_samples * (float(params["allowed_miss_pct"]) / 100), n_samples-1 )
    all_prec_quan_anno_rmNA = all_prec_quan_anno[~rm_index]
    
    rm_index = all_prec_rt_anno.loc[:, sample_names_new].isnull().values.sum(axis=1) > min( n_samples * (float(params["allowed_miss_pct"]) / 100), n_samples-1 )
    all_prec_rt_anno_rmNA = all_prec_rt_anno[~rm_index]
	# ----precur----
    
    # add Average Intensity information
    mtx0 = all_prot_quan_anno_rmNA.loc[:, sample_names_new]
    mtx0['Average Intensity'] = mtx0.mean(axis=1)
    mtx0['Average Intensity'] = mtx0['Average Intensity'].apply(np.exp2) - 1
    all_prot_quan_anno_rmNA = pd.concat([all_prot_quan_anno_rmNA.drop(columns=sample_names_new), mtx0], axis=1)
    
    # **
    all_prot_quan_anno_rmNA = reorder_dataframe(all_prot_quan_anno_rmNA,'Average Intensity','Accession.Number',1)
    
    print("\nRemoving proteins with missing values...")
    print("  # remaining proteins = " + str(all_prot_quan_anno_rmNA.shape[0]))
    all_prot_quan_anno_rmNA_csv = pd.concat([all_prot_quan_anno_rmNA.drop(columns=sample_names_new), all_prot_quan_anno_rmNA.loc[:, sample_names_new].apply(np.exp2) - 1], axis=1)
    all_prot_quan_anno_rmNA_csv.to_csv(os.path.join(interDir, "all_prot_quan_rmNA.csv"), index=False)
    # ----precur----
    mty0 = all_prec_quan_anno_rmNA.loc[:, sample_names_new]
    mty0['Average Intensity'] = mty0.mean(axis=1)
    mty0['Average Intensity'] = mty0['Average Intensity'].apply(np.exp2) - 1
    all_prec_quan_anno_rmNA = pd.concat([all_prec_quan_anno_rmNA.drop(columns=sample_names_new), mty0], axis=1)
    
    # **
    all_prec_quan_anno_rmNA = reorder_dataframe(all_prec_quan_anno_rmNA,'Average Intensity','Accession.Number',1)
    
    #print("\nRemoving precursors with missing values...")
    #print("  # remaining precursors = " + str(all_prec_quan_anno_rmNA.shape[0]))
    all_prec_quan_anno_rmNA_csv = pd.concat([all_prec_quan_anno_rmNA.drop(columns=sample_names_new), all_prec_quan_anno_rmNA.loc[:, sample_names_new].apply(np.exp2) - 1], axis=1)
    all_prec_quan_anno_rmNA_csv.to_csv(os.path.join(interDir, "all_prec_quan_rmNA.csv"), index=False)
    
    mty1 = all_prec_rt_anno_rmNA.loc[:, sample_names_new]
    mty1['Average RT'] = mty1.mean(axis=1)
    # mty1['Average RT'] = mty1['Average RT'].apply(np.exp2) - 1
    all_prec_rt_anno_rmNA = pd.concat([all_prec_rt_anno_rmNA.drop(columns=sample_names_new), mty1], axis=1)
    
    # **
    all_prec_rt_anno_rmNA = reorder_dataframe(all_prec_rt_anno_rmNA,'Average RT','Accession.Number',1)
    
    #print("\nRemoving precursors with missing values...")
    #print("  # remaining precursors = " + str(all_prec_rt_anno_rmNA.shape[0]))
    all_prec_rt_anno_rmNA_csv = all_prec_rt_anno_rmNA# pd.concat([all_prec_rt_anno_rmNA.drop(columns=sample_names_new), all_prec_rt_anno_rmNA.loc[:, sample_names_new].apply(np.exp2) - 1], axis=1)
    all_prec_rt_anno_rmNA_csv.to_csv(os.path.join(interDir, "all_prec_rt_rmNA.csv"), index=False)
	# ----precur----

    # get unique proteins
    if all_prot_quan_anno['Genes'].isnull().sum() != all_prot_quan_anno.shape[0]:
        uni_prot_quan_anno = all_prot_quan_anno[~all_prot_quan_anno['Genes'].str.contains(";", na=False)]
        uni_prot_quan_anno_csv = pd.concat([uni_prot_quan_anno.drop(columns=sample_names_new), uni_prot_quan_anno.loc[:, sample_names_new].apply(np.exp2) - 1], axis=1)
        uni_prot_quan_anno_csv.to_csv(os.path.join(pubDir, "uni_prot_quan.csv"), index=False)
        uni_prot_quan_anno_rmNA = all_prot_quan_anno_rmNA[~all_prot_quan_anno_rmNA['Genes'].str.contains(";", na=False)]
        uni_prot_quan_anno_rmNA_csv = pd.concat([uni_prot_quan_anno_rmNA.drop(columns=sample_names_new), uni_prot_quan_anno_rmNA.loc[:, sample_names_new].apply(np.exp2) - 1], axis=1)
        uni_prot_quan_anno_rmNA_csv.to_csv(os.path.join(pubDir, "uni_prot_quan_rmNA.rawinten.csv"), index=False)
        uni_prot_quan_anno_rmNA_log2 = uni_prot_quan_anno_rmNA.copy(deep=False)# uni_prot_quan_rmNA_log2
        uni_prot_quan_anno_rmNA_log2['Average Intensity'] = np.log2(uni_prot_quan_anno_rmNA_log2['Average Intensity'] + 1)
        uni_prot_quan_anno_rmNA_log2.to_csv(os.path.join(pubDir, "uni_prot_quan_rmNA.csv"), index=False)
    
    #####################################
    # Show the loading-bias information #
    #####################################
    mtx = all_prot_quan_anno_rmNA.loc[:, sample_names_new]
    new_columns = list(sample_names_new)
    new_columns.append("Protein.Names")
    mtx_PN = all_prot_quan_anno_rmNA.loc[:, new_columns]
    print("\nMedian log2(Intensity) before normalization: ")
    print(mtx.median().to_frame(name='Median of log2(Intensity)').to_string(index=True))
    
    loadBias = getLoadingBias(mtx, mtx_PN, params)
    print("\nLoading bias before normalization:")
    pd.set_option('display.max_columns', None)
    print(loadBias)
    # ----precur----
    mty = all_prec_quan_anno_rmNA.loc[:, sample_names_new]
    #new_columns = list(sample_names_new)
    #new_columns.append("Protein.Names")
    mty_PN = all_prec_quan_anno_rmNA.loc[:, new_columns]
	# ----precur----
    
    # generate tables for Rshiny
    rshiny_cols = ['Accession.Number', 'Genes', 'Protein.Description']
    all_prot_quan_anno_rmNA_rs = pd.concat([all_prot_quan_anno_rmNA.loc[:, rshiny_cols], all_prot_quan_anno_rmNA.loc[:, sample_names_new].apply(np.exp2) - 1], axis=1)
    uni_prot_quan_anno_rmNA_rs = pd.concat([uni_prot_quan_anno_rmNA.loc[:, rshiny_cols], uni_prot_quan_anno_rmNA.loc[:, sample_names_new].apply(np.exp2) - 1], axis=1)
    
    all_prot_quan_anno_rmNA_rs.to_csv(os.path.join(interDir, "all_prot_quan_rmNA.Rshiny.csv"), index=False)
    uni_prot_quan_anno_rmNA_rs.to_csv(os.path.join(pubDir, "uni_prot_quan_rmNA.Rshiny.csv"), index=False)
    # ----precur----
    rshiny_cols_prec = ['Precursor.Id', 'Genes', 'Protein.Names']
    all_prec_quan_anno_rmNA_rs = pd.concat([all_prec_quan_anno_rmNA.loc[:, rshiny_cols_prec], all_prec_quan_anno_rmNA.loc[:, sample_names_new].apply(np.exp2) - 1], axis=1)
    all_prec_quan_anno_rmNA_rs.to_csv(os.path.join(interDir, "all_prec_quan_rmNA.Rshiny.csv"), index=False)
	# ----precur----
    
    # **
    all_prot_quan_anno_rmNA.drop(columns=['Average Intensity'], inplace = True)
    # ----precur----
    all_prec_quan_anno_rmNA.drop(columns=['Average Intensity'], inplace = True)
	# ----precur----
    # Normalization
    if params["perform_normalization"] == "0":
        print("\nNormalization will not be performed!") 
        
    else: 
        if len(sample_names_new) == 1:
            print("\nOnly one sample detected. Normalization will not be performed!") 
        else:
            print("\nPerforming normalization...")
            print("  normalization method: " + params["normalization_method"])
            print("  percentage of proteins trimmed: " + params["percentage_trimmed"] + "%")
            
            #mtx_norm = global_norm(mtx, params["normalization_method"], float(params["protein_percentage"]))
            mtx_norm = normalization(mtx, mtx_PN, params)
            
            print("\nMedian log2(Intensity) after normalization: ")
            print(mtx_norm.median().to_frame(name='Median of log2(Intensity)').to_string(index=True))
            
            print("\nCalculating average intensity for each protein...")
            mtx_norm['Average Intensity'] = mtx_norm.mean(axis=1)
            
            q25 = mtx_norm['Average Intensity'].quantile(0.25) 
            q75 = mtx_norm['Average Intensity'].quantile(0.75)
            cutoff = q25 - 1.5*(q75-q25)
            print("\nRecommanded protein intensity cutoff (Q1 - 1.5 * IQR): {}, {} (log2).".format(np.exp2(cutoff), cutoff))
            
            # All
            all_prot_quan_anno_rmNA_norm = pd.concat([all_prot_quan_anno_rmNA.drop(columns=sample_names_new), mtx_norm.apply(np.exp2) - 1], axis=1)
            # **
            all_prot_quan_anno_rmNA_norm = reorder_dataframe(all_prot_quan_anno_rmNA_norm,'Average Intensity','Accession.Number',1)
            all_prot_quan_anno_rmNA_norm.to_csv(os.path.join(interDir, "all_prot_quan_rmNA_norm.csv"), index=False)
            # for Rshiny
            all_prot_quan_anno_rmNA_norm_rs = pd.concat([all_prot_quan_anno_rmNA_norm.loc[:, rshiny_cols], all_prot_quan_anno_rmNA_norm.loc[:, sample_names_new]], axis=1)
            all_prot_quan_anno_rmNA_norm_rs.to_csv(os.path.join(interDir, "all_prot_quan_rmNA_norm.Rshiny.csv"), index=False)
            # ----precur----
            mty_norm = normalization(mty, mty_PN, params)
            mty_norm['Average Intensity'] = mty_norm.mean(axis=1)
            
            # All
            all_prec_quan_anno_rmNA_norm = pd.concat([all_prec_quan_anno_rmNA.drop(columns=sample_names_new), mty_norm.apply(np.exp2) - 1], axis=1)
            # **
            all_prec_quan_anno_rmNA_norm = reorder_dataframe(all_prec_quan_anno_rmNA_norm,'Average Intensity','Accession.Number',1)
            all_prec_quan_anno_rmNA_norm.to_csv(os.path.join(interDir, "all_prec_quan_rmNA_norm.csv"), index=False)
            # for Rshiny
            all_prec_quan_anno_rmNA_norm_rs = pd.concat([all_prec_quan_anno_rmNA_norm.loc[:, rshiny_cols_prec], all_prec_quan_anno_rmNA_norm.loc[:, sample_names_new]], axis=1)
            all_prec_quan_anno_rmNA_norm_rs.to_csv(os.path.join(interDir, "all_prec_quan_rmNA_norm.Rshiny.csv"), index=False)
	        # ----precur----

            # Unique
            if all_prot_quan_anno['Genes'].isnull().sum() != all_prot_quan_anno.shape[0]:
                uni_prot_quan_anno_rmNA_norm = all_prot_quan_anno_rmNA_norm[~all_prot_quan_anno_rmNA_norm['Genes'].str.contains(";", na=False)]
                uni_prot_quan_anno_rmNA_norm.to_csv(os.path.join(pubDir, "uni_prot_quan_rmNA_norm.rawinten.csv"), index=False)
                uni_prot_quan_anno_rmNA_norm_log2 = uni_prot_quan_anno_rmNA_norm.copy(deep=False)# uni_prot_quan_rmNA_norm_log2
                uni_prot_quan_anno_rmNA_norm_log2['Average Intensity'] = np.log2(uni_prot_quan_anno_rmNA_norm_log2['Average Intensity'] + 1)
                uni_prot_quan_anno_rmNA_norm_log2_csv = pd.concat([uni_prot_quan_anno_rmNA_norm_log2.drop(columns=sample_names_new), np.log2(uni_prot_quan_anno_rmNA_norm_log2.loc[:, sample_names_new] + 1)], axis=1)
                uni_prot_quan_anno_rmNA_norm_log2_csv.to_csv(os.path.join(pubDir, "uni_prot_quan_rmNA_norm.csv"), index=False)
                # for Rshiny
                uni_prot_quan_anno_rmNA_norm_rs = pd.concat([uni_prot_quan_anno_rmNA_norm.loc[:, rshiny_cols], uni_prot_quan_anno_rmNA_norm.loc[:, sample_names_new]], axis=1)
                uni_prot_quan_anno_rmNA_norm_rs.to_csv(os.path.join(pubDir, "uni_prot_quan_rmNA_norm.Rshiny.csv"), index=False)
    
            print("\nNote: The intensity is log2 transformed in 2 files: 'uni_prot_quan_rmNA.csv' and 'uni_prot_quan_rmNA_norm.csv', others are raw intensities!")
    
    print("\nFinished!\n")






