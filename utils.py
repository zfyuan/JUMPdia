import re, sys, traceback
import numpy as np, pandas as pd
from Bio import SeqIO
import os

def getParams(paramFile):
    parameters = dict()
    with open(paramFile, 'r') as file:
        for line in file:
            #if re.search(r'^#', line) or re.search(r'^\s', line):
            if (line.strip().startswith("#") == True) or (line.strip().endswith(";") == False and "=" not in line.strip()):
                continue
            line = re.sub(r'#.*', '', line)  # Remove comments (start from '#')
            #line = re.sub(r'\s*', '', line)  # Remove all whitespaces
            line = line.strip()

            # Exception for "input_runs" parameter
            if line.endswith(";") and "input_runs" in parameters:
                parameters["input_runs"].append(line.replace(";", ""))
            else:
                key = line.split('=')[0].strip()
                val = line.split('=')[1].strip()
                if key == "input_runs":
                    parameters[key] = [val.replace(";", "")]
                else:
                    parameters[key] = val

    return parameters

# fasta_file = '/research/rgs01/home/clusterHome/yfu/projects/Core_services/diann_commandLine/proteome_fastas/human/jump_HUMAN_6_sp_tr_uc_snp_con_cu.fasta'
def extract_fasta_header(fasta_file):
    
    prot_info = pd.DataFrame(columns = ['Accession.Number','Protein.Description'])
    
    with open(fasta_file) as fasta:
        
        records = list(SeqIO.parse(fasta, 'fasta'))
        
        for record in records:
            header = record.description
            temp = header.split(" ", 1)
            temp_df = pd.DataFrame(data = {'Accession.Number': temp[0], 
                                           'Protein.Description': temp[1] if len(temp) > 1 else ''}, 
                                   index = [0])
            prot_info = pd.concat([prot_info, temp_df], ignore_index = True)
            prot_info.reset_index()
        
    return prot_info


def extract_prot_info_from_FASTAs(fasta_files, out_file):

    if len(fasta_files) == 1:
        prot_info = extract_fasta_header(fasta_files[0])
        
    elif len(fasta_files) > 1:
        prot_info = pd.DataFrame()
        for i in range(len(fasta_files)):
          temp_df = extract_fasta_header(fasta_files[i])
          prot_info = pd.concat([prot_info, temp_df], ignore_index = True)
    
    prot_info_sub = prot_info['Accession.Number'].str.split("|", expand=True)
    prot_info_sub.rename(columns={0: "Database", 1: "Protein.Ids", 2: "Protein.Names"}, inplace=True)
    
    prot_info = pd.concat([prot_info_sub, prot_info], axis=1)
    
    if out_file != 0:
        prot_info.to_csv(out_file, index = False)
    
    return prot_info


def getSubset(df_PN, params):
    # Get a subset of a dataframe to calculate loading-bias information
    df = df_PN.copy(deep=False)
    
    # 1. Filter out proteins based on the intensity level
    df["species"] = df.apply(lambda x: x["Protein.Names"].split(';')[0].split('_')[-1], axis=1)
    unique_species = set()
    for species in df["species"]:
        c_species=np.array([str(species)])
        unique_species.update(c_species)
    
    if "norm_species" in params.keys() and params["norm_species"] in unique_species:
        df = df.copy().loc[df["species"] == params["norm_species"]]
    df.drop(columns=['species','Protein.Names'], inplace = True)
    
    # calculate intensity cutoff for each sample
    q25 = df.quantile(0.25) 
    q75 = df.quantile(0.75)
    cutoff = q25 - 1.5*(q75-q25)
    
    subDf = df[(df > cutoff).prod(axis=1).astype(bool)]  

    # 2. Filter out highly variant proteins in each column 
    protMean = subDf.mean(axis=1)
    subDf = subDf.subtract(protMean, axis=0)
    pctTrimmed = float(params["percentage_trimmed"])
    n = 0
    for i in range(subDf.shape[1]):
        if n == 0:
            ind = ((subDf.iloc[:,i] > subDf.iloc[:,i].quantile(pctTrimmed / 200)) &
                   (subDf.iloc[:,i] < subDf.iloc[:,i].quantile(1 - pctTrimmed / 200)))
        else:
            ind = ind & ((subDf.iloc[:,i] > subDf.iloc[:,i].quantile(pctTrimmed / 200)) &
                         (subDf.iloc[:,i] < subDf.iloc[:,i].quantile(1 - pctTrimmed / 200)))
        n += 1

    subDf = subDf.loc[ind]
    return subDf


def getLoadingBias(df, df_PN, params):
    ###########################
    # Loading-bias evaluation #
    ###########################
    subDf = getSubset(df_PN, params)
    n = len(subDf)
    sm = 2 ** subDf.mean(axis=0)    # Sample-mean values
    msm = np.mean(sm)    # Mean of sample-mean values
    avg = sm / msm * 100
    sdVal = subDf.std(axis=0)
    #sd = ((2 ** sdVal - 1) + (1 - 2 ** (-sdVal))) / 2 * 100
    sd = ((2 ** sdVal - 2 ** (-sdVal)) / 2) * 100
    sem = sd / np.sqrt(n)
    
    loadBias = pd.concat([avg, sd, sem],axis=1)
    loadBias['#Proteins'] = n
    loadBias.rename(columns={0: "Mean[%]", 1: "SD[%]", 2: "SEM[%]"}, inplace=True)
    
    return loadBias


def normalization(df, df_PN, params):
    ################################################
    # Normalization (i.e. loading-bias correction) #
    ################################################
    res = df.copy()

    # First, get a subset for calculating normalization factors (same as loading-bias calculation)
    # Note that this subset is 1) divided by row-wise mean (i.e. PSM-wise mean) and then 2) log2-transformed
    subDf = getSubset(df_PN, params)

    # Calculate normalization factors for samples (reporters)
    if params["normalization_method"] == "mean":  # Trimmed-mean
        sm = subDf.mean(axis=0)
    elif params["normalization_method"] == "median":  # Trimmed-median
        sm = subDf.median(axis=0)
    else:
        print("Normalization method not supported!")
        
    target = np.mean(sm)
    normFactors = sm - target

    # Normalize the input dataframe, df (in log2-scale and then scale-back)
    res = res - normFactors

    return res


# Normalization
def global_norm(df, method, prot_pct):
    # df = dataframe that contains quantification data
    # method = median or mean 
    # prot_pct = percentage of non-variable proteins (chose based on CV (Coefficient of Variation)) used for normalization
    
    # calculate CV
    cv = df.apply(lambda x: np.std(x) / np.mean(x), axis=1, raw=True)
    # use proteins with smaller CV
    nonVar_idx = cv < np.quantile(cv, (prot_pct / 100))
    df_sub = df[nonVar_idx]
    # calculate correction ratio
    if method == "mean":
        norm_metric = df_sub.apply(np.mean, axis=0, raw=True)
    elif method == "median":
        norm_metric = df_sub.apply(np.median, axis=0, raw=True)
    else:
        print("Normalization method not supported!")
        
    norm_ratio = np.mean(norm_metric) / norm_metric
    # get normalized matrix
    df_norm = df.multiply(norm_ratio)
    
    return df_norm


# Context manager that copies stdout and any exceptions to a log file
class Tee(object):
    def __init__(self, filename):
        self.file = open(filename, 'w')
        self.stdout = sys.stdout

    def __enter__(self):
        sys.stdout = self

    def __exit__(self, exc_type, exc_value, tb):
        sys.stdout = self.stdout
        if exc_type is not None:
            self.file.write(traceback.format_exc())
        self.file.close()

    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)

    def flush(self):
        self.file.flush()
        self.stdout.flush()

def read_report(diann_report):
    # initial diann_report_df
    data = {"File.Name": [],"Run": [],"Protein.Group": []}
    diann_report_df = pd.DataFrame(data)
    
    # check if diann_report exists
    if not os.path.isfile(diann_report):
        return diann_report_df
    
    # check the format (*.parquet or *.tsv)
    [tpath,tname] = os.path.split(diann_report)
    [fname,ext] = os.path.splitext(tname)
    
    # read diann_report to dataframe
    if ext==".parquet":
        diann_report_df = pd.read_parquet(diann_report)
    elif ext==".tsv":
        diann_report_df = pd.read_table(diann_report)
    
    # return the dataframe
    return diann_report_df

def get_uniprot(prot_ac,protein_include_list):
    if not isinstance(prot_ac, str):
        prot_ac = str(prot_ac)
    if prot_ac.find("sp")==0 or prot_ac.find("tr")==0 or prot_ac in protein_include_list:
        prot_ac2 = "uniprot"
    else:
        prot_ac2 = "others"
    return prot_ac2

def get_KRT(prot_gn):
    if prot_gn.find("KRT")==0:
        prot_gn2 = "KRT"
    else:
        prot_gn2 = "nonKRT"
    return prot_gn2

def reorder_dataframe(df,col_pre,col_next,Large_pre_MV_pre=0):
    # all_prot_quan_anno_rmNA = reorder_dataframe(all_prot_quan_anno_rmNA,'Protein.Description','Average Intensity')
    col_list = df.columns.to_list() # Get list of columns
    if col_pre not in col_list or col_next not in col_list:
        if col_pre not in col_list:
            print('no column: {}'.format(col_pre))
        if col_next not in col_list:
            print('no column: {}'.format(col_next))
        return df
    idx_p = col_list.index(col_pre) # Find the index of 'Protein.Description'
    idx_n = col_list.index(col_next) # Find the index of 'Protein.Description'
    idx = min(idx_p,idx_n)
    if idx_p>idx_n:
        if Large_pre_MV_pre==0:
            col_list.remove(col_next)
            col_list.insert(idx_p, col_next)
        else:
            col_list.remove(col_pre)
            col_list.insert(idx, col_pre)
    else:
        col_list.remove(col_next) # Move 'Average Intensity' to the position right after 'Protein.Description'
        col_list.insert(idx + 1, col_next)
    df = df[col_list] # Reorder the DataFrame
    return df
