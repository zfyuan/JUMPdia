#!/usr/bin/python

import os
import fnmatch
import string
import pickle
import shutil
import datetime

def Get_Params(paramFile):
	# get params
	dict1 = {}
	with open(paramFile,"r") as f:
		for line in f:
			if (line.startswith("#") == False) and ("=" in line.strip()):
	
				te_line = line.strip().split("=")
				key = te_line[0].strip()

				if "#" in te_line[1]:
					value = te_line[1].split("#")[0].strip()
				else:
					value = te_line[1].strip()
				dict1[key]= value
	return dict1

def Check_Path(output_path):
	# check and create path
	if False==os.path.isdir(output_path):
		os.mkdir(output_path)
		if False==os.path.isdir(output_path):
			raise(Exception,output_path)
	return

def Check_OneFormat(SourcePath,extend):
	# get the full filename and only filename of the extended
	matched_fullfiles = []
	matched_files = []
	dir_list = os.listdir(SourcePath)
	for d in dir_list:
		if os.path.isfile(os.path.join(SourcePath,d)):
			if fnmatch.fnmatch(d.lower(),extend.lower()):
				matched_fullfiles.append(os.path.join(SourcePath,d))
				matched_files.append(d)
	return matched_fullfiles,matched_files

def Check_d_Folder(SourcePath):
	# get the full filename and only filename of the extended
	extend = '*.d'
	matched_fullfiles = []
	matched_files = []
	dir_list = os.listdir(SourcePath)
	for d in dir_list:
		if os.path.isdir(os.path.join(SourcePath,d)):
			if fnmatch.fnmatch(d.lower(),extend.lower()):
				matched_fullfiles.append(os.path.join(SourcePath,d))
				matched_files.append(d)
	return matched_fullfiles,matched_files

def Get_mzXMLlist_without_ext(update_path):
	# get mzXMLlist
	# mzXMLlist_fullfile = os.path.join(update_path,'c_mzXMLlist.dat')
	# if True==os.path.isfile(mzXMLlist_fullfile):
		# mzXMLlist = pickle.load(open(mzXMLlist_fullfile,'rb'))
		# return mzXMLlist
	
	[all_fullfiles,mzXML_files] = Check_OneFormat(update_path,'*.mzXML')
	mzXML_files.sort()
	mzXMLlist = []
	for i in range(0,len(mzXML_files)):
		[fname,ext] = os.path.splitext(mzXML_files[i])
		mzXMLlist.append(fname)
	# pickle.dump(mzXMLlist,open(mzXMLlist_fullfile,'wb'))
	return mzXMLlist

def Get_mzMLlist_without_ext(update_path):
	# get mzMLlist
	# mzMLlist_fullfile = os.path.join(update_path,'c_mzMLlist.dat')
	# if True==os.path.isfile(mzMLlist_fullfile):
		# mzMLlist = pickle.load(open(mzMLlist_fullfile,'rb'))
		# return mzMLlist
	
	[all_fullfiles,mzML_files] = Check_OneFormat(update_path,'*.mzML')
	mzML_files.sort()
	mzMLlist = []
	for i in range(0,len(mzML_files)):
		[fname,ext] = os.path.splitext(mzML_files[i])
		mzMLlist.append(fname)
	# pickle.dump(mzMLlist,open(mzMLlist_fullfile,'wb'))
	return mzMLlist

def Get_rawlist_without_ext(update_path):
	# get rawlist
	# rawlist_fullfile = os.path.join(update_path,'c_rawlist.dat')
	# if True==os.path.isfile(rawlist_fullfile):
		# rawlist = pickle.load(open(rawlist_fullfile,'rb'))
		# return rawlist
	
	[all_fullfiles,raw_files] = Check_OneFormat(update_path,'*.raw')
	raw_files.sort()
	rawlist = []
	for i in range(0,len(raw_files)):
		[fname,ext] = os.path.splitext(raw_files[i])
		rawlist.append(fname)
	# pickle.dump(rawlist,open(rawlist_fullfile,'wb'))
	return rawlist

def Get_dFolderlist(update_path):
	# get dFolderlist
	# dFolderlist_fullfile = os.path.join(update_path,'c_dFolderlist.dat')
	# if True==os.path.isfile(dFolderlist_fullfile):
		# dFolderlist = pickle.load(open(dFolderlist_fullfile,'rb'))
		# return dFolderlist
	
	# [all_fullfiles,dFolder_files] = Check_d_Folder(update_path)
	# dFolder_files.sort()
	# dFolderlist = []
	# for i in range(0,len(dFolder_files)):
		# fname = dFolder_files[i]#[fname,ext] = os.path.splitext(dFolder_files[i])
		# dFolderlist.append(fname)
	[all_fullfiles,dFolderlist] = Check_d_Folder(update_path)
	dFolderlist.sort()
	# pickle.dump(dFolderlist,open(dFolderlist_fullfile,'wb'))
	return dFolderlist

def GetUpdateDate():
	# get update date of today
	update1 = datetime.date.today() # today
	
	update2 = '%s' % update1 # obj to string, e.g. 2020-04-22
	update = '%s%s%s' % (update2[0:4],update2[5:7],update2[8:10]) # string, e.g. 20200420
	return update

def GetFastaPathName(fasta_file):
	# get fasta path and name
	fastas = fasta_file.split(";")
	fasta_name = []
	if len(fastas)==1:
		[fasta_path,tname] = os.path.split(fastas[0])
		#[c_fasta_name,ext] = os.path.splitext(tname)
		fasta_name.append(tname)
	else:
		fastas.sort()
		for fasta in fastas:
			[fasta_path,tname] = os.path.split(fasta)
			#[c_fasta_name,ext] = os.path.splitext(tname)
			fasta_name.append(tname)
	
	[fasta_path,tname] = os.path.split(fastas[0])
	
	return fasta_path,fasta_name

def GetFastaShortName(fasta_file,params):
	# get fasta short name
	cmd_prt1 = "" # for modifications
	enzyme_short = {'Trypsin_P': 'TryP',
				'Trypsin': 'Tryp',
				'Lys_C': 'LysC',
				'Chymotrypsin': 'Chym',
				'AspN': 'AspN',
				'GluC': 'GluC',}
	cmd_prt1 = cmd_prt1 + enzyme_short[params['protease']] + "_mc" + params['n_missed_cleavages']
	if params['C_carbamidomethylation'] == '1':
		cmd_prt1 = cmd_prt1 + "_c57"
	elif params['C_carbamidomethylation'] == '0':
		cmd_prt1 = cmd_prt1 + "_c0"
	if params['Phosphorylation'] == '1':
		cmd_prt1 = cmd_prt1 + "_pho"
	if params['N_term_acetylation'] == '1' or params['K_GG_addcut'] == '1' or '--var-mod' in params['additional_options'] or '--fixed-mod' in params['additional_options']:
		if '--var-mod SILAC' in params['additional_options']:
			mod3 = 'SILAC1var'
		elif '--fixed-mod SILAC' in params['additional_options']:
			mod3 = 'SILAC2fix'
		elif '--var-mod' in params['additional_options']:
			mod3 = 'ADD1'
		else:
			mod3 = 'ADD0'
		cmd_prt1 = cmd_prt1 + "_mod" + params['N_term_acetylation'] + params['K_GG_addcut'] + mod3
	
	fastas = fasta_file.split(";")
	if len(fastas)==1:
		[fasta_path,tname] = os.path.split(fastas[0])
		[c_fasta_name,ext] = os.path.splitext(tname)
		fasta_short_name = c_fasta_name.split('_')[0].lower() + "_" + cmd_prt1
	else:
		fastas.sort()
		mix_fasta_name = ""
		for fasta in fastas:
			[fasta_path,tname] = os.path.split(fasta)
			[c_fasta_name,ext] = os.path.splitext(tname)
			mix_fasta_name = mix_fasta_name + c_fasta_name[0].lower()
		fasta_short_name = mix_fasta_name + "_" + cmd_prt1
	
	return fasta_short_name

def generate_diann_LibGen_params(cur_path,data_path,fasta_file,dFolderlist,mzMLlist,rawlist,update_date,params):
	# generate diann_LibGen params
	
	c_diann_LibGen_fullfile = os.path.join(cur_path,'diann_1_LibGen.params')
	#if True==os.path.isfile(c_diann_LibGen_fullfile):
	#	continue
	
	# fasta_short_name
	fasta_short_name = GetFastaShortName(fasta_file,params)
	fastas = fasta_file.split(";")
	
	# open, write, close
	file1 = open(c_diann_LibGen_fullfile,'w')
	
	file1.write('# DIA-NN in silico spectral library generation parameters (Version: 0.1.0, Date: 10/08/2025)\n')
	file1.write('\n')
	file1.write('# FASTA files (use ";" to combine multiple FASTA files to one row for "fasta_file")\n')
	if len(fastas)==1:
		file1.write('fasta_file = {}\n'.format(fastas[0]))
	else:
		file1.write('fasta_file = {}'.format(fastas[0]))
		for i in range(1,len(fastas)):
			file1.write(';{}\n'.format(fastas[i]))
	file1.write('\n')
	file1.write('# Output library name\n')
	file1.write('library_name = {}_lib                    # the library will be saved in .speclib format\n'.format(fasta_short_name))
	file1.write('\n')
	file1.write('# Digestion rule\n')
	file1.write('protease = {}                                  # default: Trypsin_P (involve cuts at K* and R*)\n'.format(params['protease']))
	file1.write('                                                      # other options: Trypsin (involve cuts at K* and R*, but excluding cuts at *P)\n')
	file1.write('                                                      #                Lys_C (involve cuts at K*)\n')
	file1.write('                                                      #                Chymotrypsin (involve cuts at F*, Y*, W*, M* and L*, but excluding cuts at *P)\n')
	file1.write('                                                      #                AspN (involve cuts at *D)\n')
	file1.write('                                                      #                GluC (involve cuts at E*)\n')
	file1.write('n_missed_cleavages = {}                                # number of missed cleavages\n'.format(params['n_missed_cleavages']))
	file1.write('min_pep_length = {}                                    # minimum peptide length\n'.format(params['min_pep_length']))
	file1.write('max_pep_length = {}                                   # maximum peptide length\n'.format(params['max_pep_length']))
	file1.write('\n')
	file1.write('# Fixed modifications\n')
	file1.write('C_carbamidomethylation = {}                            # enable cystine carbamidomethylation as a fixed modification, 0 = off; 1 = on\n'.format(params['C_carbamidomethylation']))
	file1.write('\n')
	file1.write('# Variable modifications\n')
	file1.write('M_oxidation = {}                                       # enable methionine oxidation as a variable modification, 0 = off; 1= on\n'.format(params['M_oxidation']))
	file1.write('N_term_M_excision = {}                                 # enable protein N-terminal methionine excision as a variable modification, 0 = off; 1 = on\n'.format(params['N_term_M_excision']))
	file1.write('max_n_variable_mod = {}                                # maximum number of variable modifications\n'.format(params['max_n_variable_mod']))
	file1.write('\n')
	file1.write('# Post-translational modifications (PTMs) \n')
	file1.write('N_term_acetylation = {}                                # enable protein N-terminal acetylation as a variable modification / score sites, 0 = off; 1 = on\n'.format(params['N_term_acetylation']))
	file1.write('Phosphorylation = {}                                   # enable S, T, Y phosphorylation as a variable modification / score sites, 0 = off; 1 = on\n'.format(params['Phosphorylation']))
	file1.write('K_GG_addcut = {}                                       # enable -GG adduct on lysines (left after tryptic digest of an attached ubiquitin) as a variable modification / score sites, 0 = off; 1 = on\n'.format(params['K_GG_addcut']))
	file1.write('\n')
	file1.write('# Precursor ion properties\n')
	file1.write('min_prec_charge = {}                                   # minimum precursor charge\n'.format(params['min_prec_charge']))
	file1.write('max_prec_charge = {}                                   # maximum precursor charge\n'.format(params['max_prec_charge']))
	file1.write('min_prec_mz = {}                                     # minimum precursor m/z\n'.format(params['min_prec_mz']))
	file1.write('max_prec_mz = {}                                    # maximum precursor m/z\n'.format(params['max_prec_mz']))
	file1.write('\n')
	file1.write('# Fragment ion properties\n')
	file1.write('min_frag_mz = {}                                     # minimum fragment ion m/z\n'.format(params['min_frag_mz']))
	file1.write('max_frag_mz = {}                                    # maximum fragment ion m/z\n'.format(params['max_frag_mz']))
	file1.write('\n')
	file1.write('# Other parameters\n')
	file1.write('n_CPU_threads = {}                                    # number of CPU threads used for running\n'.format(params['n_CPU_threads']))
	file1.write('log_level = {}                                         # specify how detailed the information reported in the log file, use numbers between 0~5\n'.format(params['log_level']))
	file1.write('additional_options = {}                             # null (default); others, e.g. --var-mod Dimethyl,28.0313,KR\n'.format(params['additional_options']))
	
	file1.close()
	return

def generate_diann_Search_params(cur_path,data_path,fasta_file,dFolderlist,mzMLlist,rawlist,update_date,params):
	# generate diann_Search params
	
	# data_path_ln (softlink of data_path)
	data_path_ln = os.path.join( os.path.split(cur_path)[0],'data')
	
	c_diann_Search_fullfile = os.path.join(cur_path,'diann_2_Search.params')
	#if True==os.path.isfile(c_diann_Search_fullfile):
	#	continue
	
	# fasta_path,fasta_name
	[fasta_path,fasta_name] = GetFastaPathName(fasta_file)
	
	# fasta_short_name
	fasta_short_name = GetFastaShortName(fasta_file,params)
	fastas = fasta_file.split(";")
	
	# open, write, close
	file1 = open(c_diann_Search_fullfile,'w')
	
	file1.write('# DIA-NN search and quantification parameters (Version: 0.1.0, Date: 10/08/2025)\n')
	file1.write('\n')
	file1.write('# Input: mzML or raw files (for Thermo instruments) or diaPASEF .d files (for timsTOF instruments)\n')
	file1.write('# for multiple runs, end each file with ";"\n')
	if len(dFolderlist)>0:
		file1.write('input_runs = {};\n'.format(os.path.join(data_path_ln,dFolderlist[0])))
		if len(dFolderlist)>1:
			for i in range(1,len(dFolderlist)):
				file1.write('{};\n'.format(os.path.join(data_path_ln,dFolderlist[i])))
	elif len(mzMLlist)>0:
		file1.write('input_runs = {}.mzML;\n'.format(os.path.join(data_path_ln,mzMLlist[0])))
		if len(mzMLlist)>1:
			for i in range(1,len(mzMLlist)):
				file1.write('{}.mzML;\n'.format(os.path.join(data_path_ln,mzMLlist[i])))
	elif len(rawlist)>0:
		file1.write('input_runs = {}.raw;\n'.format(os.path.join(data_path_ln,rawlist[0])))
		if len(rawlist)>1:
			for i in range(1,len(rawlist)):
				file1.write('{}.raw;\n'.format(os.path.join(data_path_ln,rawlist[i])))
	else:
		file1.write('input_runs = ;\n')
	file1.write('\n')
	file1.write('# Spectral library (recommend using the .speclib format file)\n')
	file1.write('SpecLib = {}\n'.format(os.path.join(fasta_path,fasta_short_name+'_lib',fasta_short_name+'_lib.predicted.speclib')))
	file1.write('\n')
	file1.write('# Fasta files (use ";" to combine multiple FASTA files to one row for "fasta_file")\n')
	if len(fastas)==1:
		file1.write('fasta_file = {}\n'.format(fastas[0]))
	else:
		file1.write('fasta_file = {}'.format(fastas[0]))
		for i in range(1,len(fastas)):
			file1.write(';{}\n'.format(fastas[i]))
	file1.write('re_annotate = 0                                       # 0 = false (default), 1 = true\n')
	file1.write('\n')
	file1.write('# Use existing .quant files\n')
	file1.write('use_quant = 0                                         # 0 = false (default), 1 = true\n')
	file1.write('quantities_matrices = 1                               # 0 = false, 1 = true (default)\n')
	file1.write('XICs = 1                                              # 0 = false, 1 = true (default)\n')
	file1.write('\n')
	file1.write('\n')
	file1.write('# Output directory name\n')
	file1.write('outputDir = search_results_{}\n'.format(update_date))
	file1.write('\n')
	file1.write('# Identification related parameters\n')
	file1.write('ms2_tolerance = {}                                     # MS2 mass accuracy (ppm), use 0 for automatic inference\n'.format(params['ms2_tolerance']))
	file1.write('ms1_tolerance = {}                                     # MS1 mass accuracy (ppm), use 0 for automatic inference\n'.format(params['ms1_tolerance']))
	file1.write('n_scans_per_peak = {}                                  # number of scans per peak, use 0 for automatic inference\n'.format(params['n_scans_per_peak']))
	#file1.write('use_isotopologues = 1                                 # extract precursor isotopologues, 0 = off; 1 = on\n')
	#file1.write('noise_precursor_ID_removal = 1                        # corresponding to "No shared spectra" in the GUI, 0 = off; 1 = on\n')
	file1.write('unrelated_runs = {}                                    # Different runs will be treated as unrelated, 0 = false (default); 1 = true\n'.format(params['unrelated_runs']))
	file1.write('fdr_prec = {}                                       # FDR at precursor level\n'.format(params['fdr_prec']))
	file1.write('scoring = {}                                           # 1 = Generic (default); 2 = Peptidoforms; 3 = Proteoforms\n'.format(params['scoring']))
	file1.write('machine_learning = {}                                  # 1 = Linear classifier; 2 = NNs (fast); 3 = NNs (cross-validated) (default)\n'.format(params['machine_learning']))
	file1.write('\n')
	file1.write('# Fixed modifications\n')
	file1.write('C_carbamidomethylation = {}                            # enable cystine carbamidomethylation as a fixed modification, 0 = off; 1 = on\n'.format(params['C_carbamidomethylation']))
	file1.write('\n')
	file1.write('# Variable modifications\n')
	file1.write('M_oxidation = {}                                       # enable methionine oxidation as a variable modification, 0 = off; 1= on\n'.format(params['M_oxidation']))
	file1.write('N_term_M_excision = {}                                 # enable protein N-terminal methionine excision as a variable modification, 0 = off; 1 = on\n'.format(params['N_term_M_excision']))
	file1.write('max_n_variable_mod = {}                                # maximum number of variable modifications\n'.format(params['max_n_variable_mod']))
	file1.write('\n')
	file1.write('# Post-translational modifications (PTMs) \n')
	file1.write('N_term_acetylation = {}                                # enable protein N-terminal acetylation as a variable modification / score sites, 0 = off; 1 = on\n'.format(params['N_term_acetylation']))
	file1.write('Phosphorylation = {}                                   # enable S, T, Y phosphorylation as a variable modification / score sites, 0 = off; 1 = on\n'.format(params['Phosphorylation']))
	file1.write('K_GG_addcut = {}                                       # enable -GG adduct on lysines (left after tryptic digest of an attached ubiquitin) as a variable modification / score sites, 0 = off; 1 = on\n'.format(params['K_GG_addcut']))
	file1.write('\n')
	file1.write('# Match between runs (MBR)\n')
	file1.write('MBR = {}                                               # create a spectral library from DIA data and reanalyze using this library, 0 = off; 1 = on (default)\n'.format(params['MBR']))
	file1.write('\n')
	file1.write('# Quantification related parameters\n')
	#file1.write('frag_noise_removal = 1                                # interference subtraction from fragment ion elution peaks, 0 = off (high precision); 1 = on (high accuracy)\n')
	#file1.write('robust_LC = 0                                         # 0 = Any LC; 1 = Robust LC, "Robust LC" instructs DIA-NN to only integrate peaks with close apex RT\n')
	#file1.write('use_peak_hight = 0                                    # use the apex height of the peak for quantification instead of area under the peak, 0 = off; 1 = on\n')
	file1.write('quantification_strategy = {}                           # 1 = Legacy (direct); 2 = QuantUMS (high accuracy); 3 = QuantUMS (high precision) (default)\n'.format(params['quantification_strategy']))
	file1.write('cross_run_normalization = {}                           # 1 = Global; 2 = RT-dependent (default); 3 = RT & signal-dep. (experimental); 4 = off\n'.format(params['cross_run_normalization']))
	file1.write('\n')
	file1.write('# Protein inference mode\n')
	file1.write('prot_inference = {}                                    # 0 = off; 1 = on (default), if turned on, one protein will be assigned to only one protein group\n'.format(params['prot_inference']))
	file1.write('proteotypicity = {}                                    # 1 = Isofrom IDs; 2 = Protein names (from FASTA); 3 = Genes (Species-specific); 4 = Genes (default)\n'.format(params['proteotypicity']))
	file1.write('\n')
	file1.write('# Other parameters\n')
	#file1.write('report_LibInfo = 0                                    # add extra library information on the precursor and its fragments to the main output report, 0 = off; 1 = on\n')
	file1.write('library_generation = {}                                # 1 = IDs profiling; 2 = IDs, RT & IM profiling (default); 3 = Smart profiling; 4 = Full profiling\n'.format(params['library_generation']))
	file1.write('speed_RAM_usage = {}                                   # 1 = Optimal results (default); 2 = Low RAM usage; 3 = Low RAM & high speed; 4 = Ultra-fast\n'.format(params['speed_RAM_usage']))
	file1.write('n_CPU_threads = {}                                    # number of CPU threads used for running\n'.format(params['n_CPU_threads']))
	file1.write('log_level = {}                                         # specify how detailed the information reported in the log file, use numbers between 0~5\n'.format(params['log_level']))
	file1.write('additional_options = {}                             # null (default); others, e.g. --var-mod Dimethyl,28.0313,KR\n'.format(params['additional_options']))
	
	file1.close()
	return

def generate_diann_iqQuan_params(cur_path,data_path,fasta_file,dFolderlist,mzMLlist,rawlist,update_date,params):
	# generate diann_iqQuan params
	
	# data_path_ln (softlink of data_path)
	data_path_ln = os.path.join( os.path.split(cur_path)[0],'data')
	
	c_diann_iqQuan_fullfile = os.path.join(cur_path,'diann_3_iqQuan.params')
	#if True==os.path.isfile(c_diann_iqQuan_fullfile):
	#	continue
	
	# fasta_path,fasta_name
	[fasta_path,fasta_name] = GetFastaPathName(fasta_file)
	
	# fasta_short_name
	fasta_short_name = GetFastaShortName(fasta_file,params)
	
	# open, write, close
	file1 = open(c_diann_iqQuan_fullfile,'w')
	
	file1.write('# Protein-level quantification parameters (Version: 0.1.0, Date: 10/08/2025)\n')
	file1.write('\n')
	file1.write('# DIA-NN report file (report.parquet)\n')
	file1.write('diann_report = {}\n'.format(os.path.join(data_path_ln,'search_results_'+update_date,'report.parquet')))
	file1.write('\n')
	file1.write('# Protein annotation file (This can be generated by paring the FASTA file)\n')
	file1.write('prot_anno_file = {}\n'.format(os.path.join(fasta_path,fasta_short_name+'_prot_anno',fasta_short_name+'_prot_anno.csv')))
	file1.write('\n')
	file1.write('# Output directory\n')
	file1.write('outputDir = quan_results_{}\n'.format(update_date))
	file1.write('\n')
	file1.write('# Quantification method (Important!)\n')
	file1.write('QuanMethod = {}                                        # 0 = DIANN native (use PG.TopN)\n'.format(params['QuanMethod']))
	file1.write('                                                      # 1 = MaxLFQ (use PG.MaxLFQ)\n')
	file1.write('\n')
	file1.write('# Q value cutoffs\n')
	file1.write('Lib.Q.Value = {}                                    # Please use the same column names from DIA-NN report\n'.format(params['Lib.Q.Value']))
	file1.write('Lib.PG.Q.Value = {}\n'.format(params['Lib.PG.Q.Value']))
	file1.write('\n')
	file1.write('# Filter\n')
	file1.write('allowed_miss_pct = {}                                # filter out proteins based on missing value percentages\n'.format(params['allowed_miss_pct']))
	file1.write('\n')
	file1.write('# Normalization\n')
	file1.write('perform_normalization = {}                             # 0 = off, 1 = on\n'.format(params['perform_normalization']))
	file1.write('normalization_method = {}                         # use median or mean intensity for cross-run normalization\n'.format(params['normalization_method']))
	file1.write('percentage_trimmed = {}                               # percentage of most variable intensities to be trimmed for the normalization\n'.format(params['percentage_trimmed']))
	file1.write('norm_species = {}                                    # if it is on, normalized to : ALL (default), 1 species for Benchmark: HUMAN; MOUSE; RAT; YEAST; SCHPO; ECOLI; DROME; ARATH; CAEEL; CHICK; DANRE\n'.format(params['norm_species']))
	
	file1.close()
	return

def generate_proteoDA_files(cur_path,data_path,fasta_file,dFolderlist,mzMLlist,rawlist,update_date,params):
	# generate proteoDA files
	
	proteoDA_params_R = '/hpcf/authorized_apps/proteomics_apps/proteoDA2/proteoDA_params.R'
	shutil.copy(proteoDA_params_R, cur_path)
	
	# data_path_ln (softlink of data_path)
	#proteoDA2_path = os.path.join(os.path.split(cur_path)[0],'proteoDA2')
	
	c_contrasts_fullfile = os.path.join(cur_path,'proteoDA_contrasts.csv')
	c_metadata_fullfile = os.path.join(cur_path,'proteoDA_Sample_metadata.csv')
	#if True==os.path.isfile(c_diann_iqQuan_fullfile):
	#	continue
	
	# open, write, close
	file1 = open(c_contrasts_fullfile,'w')
	file1.write('group1_vs_group2=group1-group2\n')
	file1.close()
	
	file2 = open(c_metadata_fullfile,'w')
	file2.write('file,sample,group\n')
	if len(dFolderlist)>0:
		if len(dFolderlist)==1:
			file2.write('{},{},{}\n'.format(dFolderlist[0].replace('-','_').split('.d')[0],dFolderlist[0].replace('-','_').split('.d')[0],'group1'))
		else:
			half = int(len(dFolderlist)/2)
			for i in range(0,half):
				file2.write('{},{},{}\n'.format(dFolderlist[i].replace('-','_').split('.d')[0],dFolderlist[i].replace('-','_').split('.d')[0],'group1'))
			for i in range(half,len(dFolderlist)):
				file2.write('{},{},{}\n'.format(dFolderlist[i].replace('-','_').split('.d')[0],dFolderlist[i].replace('-','_').split('.d')[0],'group2'))
	elif len(mzMLlist)>0:
		if len(mzMLlist)==1:
			file2.write('{},{},{}\n'.format(mzMLlist[0].replace('-','_').split('.d')[0],mzMLlist[0].replace('-','_').split('.d')[0],'group1'))
		else:
			half = int(len(mzMLlist)/2)
			for i in range(0,half):
				file2.write('{},{},{}\n'.format(mzMLlist[i].replace('-','_').split('.d')[0],mzMLlist[i].replace('-','_').split('.d')[0],'group1'))
			for i in range(half,len(mzMLlist)):
				file2.write('{},{},{}\n'.format(mzMLlist[i].replace('-','_').split('.d')[0],mzMLlist[i].replace('-','_').split('.d')[0],'group2'))
	elif len(rawlist)>0:
		if len(rawlist)==1:
			file2.write('{},{},{}\n'.format(rawlist[0].replace('-','_').split('.d')[0],rawlist[0].replace('-','_').split('.d')[0],'group1'))
		else:
			half = int(len(rawlist)/2)
			for i in range(0,half):
				file2.write('{},{},{}\n'.format(rawlist[i].replace('-','_').split('.d')[0],rawlist[i].replace('-','_').split('.d')[0],'group1'))
			for i in range(half,len(rawlist)):
				file2.write('{},{},{}\n'.format(rawlist[i].replace('-','_').split('.d')[0],rawlist[i].replace('-','_').split('.d')[0],'group2'))
	else:
		file2.write(',,\n')
	file2.close()
	
	return

def generate_workflow(cur_path,data_path,fasta_file,dFolderlist,mzMLlist,rawlist,update_date,params):
	# generate jump -q bash
	
	# data_path_ln (softlink of data_path)
	data_path_ln = os.path.join( os.path.split(cur_path)[0],'data')
	
	proteoDA2_path = os.path.join(os.path.split(cur_path)[0],'proteoDA2')
	
	c_workflow_sh_fullfile = os.path.join(cur_path,'workflow.sh')
	#if True==os.path.isfile(c_workflow_sh_fullfile):
	#	continue
	
	# fasta_path,fasta_name
	[fasta_path,fasta_name] = GetFastaPathName(fasta_file)
	
	# fasta_short_name
	fasta_short_name = GetFastaShortName(fasta_file,params)
	
	# predicted_speclib_fullfile, prot_anno_fullfile
	predicted_speclib_fullfile = os.path.join(fasta_path,fasta_short_name+'_lib',fasta_short_name+'_lib'+'.predicted.speclib')
	prot_anno_fullfile = os.path.join(fasta_path,fasta_short_name+'_prot_anno',fasta_short_name+'_prot_anno'+'.csv')
	
	# open, write, close
	file1 = open(c_workflow_sh_fullfile,'w')
	
	file1.write('#!/bin/bash\n')
	file1.write('module purge\n')
	# file1.write('# load dotnet to run raw files\n')# load dotnet
	# # file1.write('rm -rf ~/.config/Spectronaut\n')# clear some metadata for Spectronaut
	# # file1.write('module load spectronaut/19.5\n')# Spectronaut version '19.5' loaded; .NET version '8.0.404' loaded
	# file1.write('module load dotnet/8.0.407\n')# .NET version '8.0.407' loaded
	file1.write('# activate the ‘diann_iq_cmd’ Conda environment\n')
	file1.write('module load gcc/13.1.0\n')
	file1.write('module load conda3/202402\n')
	file1.write('conda activate /hpcf/authorized_apps/proteomics_apps/conda/proteoda2\n')
	# file1.write('module load miniconda\n')
	# file1.write('source /hpcf/authorized_apps/proteomics_apps/DIANN_CMD/diann_iq_cmd/bin/activate\n')
	file1.write('\n')
	file1.write('# make softlink of the data_path\n')
	file1.write('mkdir -p {}\n'.format(data_path_ln))
	if len(dFolderlist)>0:
		file1.write('ln -s {} {}\n'.format(os.path.join(data_path,'*.d'),data_path_ln))
	elif len(mzMLlist)>0:
		file1.write('ln -s {} {}\n'.format(os.path.join(data_path,'*.mzML'),data_path_ln))
	elif len(rawlist)>0:
		#file1.write('ln -s {} {}\n'.format(os.path.join(data_path,'*.raw'),data_path_ln))
		file1.write('cp {} {}\n'.format(os.path.join(data_path,'*.raw'),data_path_ln))
	else:
		file1.write('# keep using data_path without softlink\n')
		data_path_ln = data_path
	file1.write('\n')
	file1.write('# copy parameter files\n')
	if os.path.isfile(predicted_speclib_fullfile)==False or os.path.isfile(prot_anno_fullfile)==False:
		file1.write('cp {} {}\n'.format(os.path.join(cur_path,'diann_1_LibGen.params'),fasta_path))
	else:
		file1.write('#cp {} {}\n'.format(os.path.join(cur_path,'diann_1_LibGen.params'),fasta_path))
	file1.write('cp {} {}\n'.format(os.path.join(cur_path,'diann_2_Search.params'),data_path_ln))
	file1.write('cp {} {}\n'.format(os.path.join(cur_path,'diann_3_iqQuan.params'),data_path_ln))
	file1.write('\n')
	if os.path.isfile(predicted_speclib_fullfile)==False or os.path.isfile(prot_anno_fullfile)==False:
		file1.write('# predict the library\n')
		file1.write('cd {}\n'.format(fasta_path))
		file1.write('/hpcf/authorized_apps/proteomics_apps/conda/proteoda2/bin/python /hpcf/authorized_apps/proteomics_apps/DIANN_CMD/diann_iq/diann2.2.0/20251118/diann_LibGen.py diann_1_LibGen.params\n')
		file1.write('/hpcf/authorized_apps/proteomics_apps/conda/proteoda2/bin/python /hpcf/authorized_apps/proteomics_apps/DIANN_CMD/diann_iq/diann2.2.0/20251118/parse_FASTA.py {} {}\n'.format(fasta_short_name+'_prot_anno',' '.join(fasta_name)))
		file1.write('\n')
	else:
		file1.write('# predict the library (skipped as they already exist)\n')
		file1.write('#cd {}\n'.format(fasta_path))
		file1.write('#/hpcf/authorized_apps/proteomics_apps/conda/proteoda2/bin/python /hpcf/authorized_apps/proteomics_apps/DIANN_CMD/diann_iq/diann2.2.0/20251118/diann_LibGen.py diann_1_LibGen.params\n')
		file1.write('#/hpcf/authorized_apps/proteomics_apps/conda/proteoda2/bin/python /hpcf/authorized_apps/proteomics_apps/DIANN_CMD/diann_iq/diann2.2.0/20251118/parse_FASTA.py {} {}\n'.format(fasta_short_name+'_prot_anno',' '.join(fasta_name)))
		file1.write('\n')
	file1.write('# search and quantify precursors, quantify proteins\n')
	file1.write('cd {}\n'.format(data_path_ln))
	file1.write('/hpcf/authorized_apps/proteomics_apps/conda/proteoda2/bin/python /hpcf/authorized_apps/proteomics_apps/DIANN_CMD/diann_iq/diann2.2.0/20251118/diann_Search.py diann_2_Search.params\n')
	file1.write('/hpcf/authorized_apps/proteomics_apps/conda/proteoda2/bin/python /hpcf/authorized_apps/proteomics_apps/DIANN_CMD/diann_iq/diann2.2.0/20251118/diann_iqQuan.py diann_3_iqQuan.params\n')
	file1.write('\n')
	
	if params['proteoDA2'] == '1':
		file1.write('# === Custom proteoDA post-processing ===\n')# proteoDA
		file1.write('# !!!generate and edit: CMD/proteoDA_params.R, CMD/proteoDA_contrasts.csv, CMD/proteoDA_Sample_metadata.csv\n')# edit proteoDA files
		generate_proteoDA_files(cur_path,data_path,fasta_file,dFolderlist,mzMLlist,rawlist,update_date,params)
		file1.write('\n')
		
		file1.write('# copy proteoDA files: from "CMD" to "proteoDA2" and "proteoDA2/data"\n')# copy proteoDA files
		file1.write('cp -R /hpcf/authorized_apps/proteomics_apps/proteoDA2/ {}\n'.format( os.path.split(cur_path)[0] ))
		file1.write('cp {} {}\n'.format( os.path.join(cur_path,'proteoDA_params.R'),proteoDA2_path ))
		file1.write('cp {} {}\n'.format( os.path.join(cur_path,'proteoDA_contrasts.csv'),os.path.join(proteoDA2_path,'data','contrasts.csv') ))
		file1.write('cp {} {}\n'.format( os.path.join(cur_path,'proteoDA_Sample_metadata.csv'),os.path.join(proteoDA2_path,'data','Sample_metadata.csv') ))
		file1.write('\n')
		
		file1.write('cd {}\n'.format(data_path_ln))
		file1.write('latest_quan_dir=$(find \"$(pwd)\" -maxdepth 1 -type d -name \"quan_*\" | sort | tail -n 1)\n')
		file1.write('echo \"[INFO] Using DIANN output from: $latest_quan_dir\"\n')
		file1.write('diann_output=\"$latest_quan_dir/publication/uni_prot_quan_rmNA_norm.csv\"\n')
		file1.write('echo \"[INFO] Running R script with input file: $diann_output\"\n')
		file1.write('cd {}\n'.format(proteoDA2_path))
		file1.write('Rscript proteoDA_diann_per_contrasts_ctrlProteins.R \"$diann_output\" > >(tee proteoDA_stdout.log) 2> >(tee proteoDA_stderr.log >&2)\n')
		file1.write('wait\n')
		file1.write('Rscript GenerateReports.R\n')
		file1.write('wait\n')
		file1.write('\n')
	
	file1.write('# deactivate the ‘diann_iq_cmd’ Conda environment after the task is complete\n')
	file1.write('conda deactivate\n')
	# file1.write('source /hpcf/authorized_apps/proteomics_apps/DIANN_CMD/diann_iq_cmd/bin/deactivate\n')
	
	file1.close()
	return
