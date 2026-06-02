#!/usr/bin/python

import os
import sys
import time
import platform

import DIANNParamsLib

def main():

# start #
	# start_time = time.time()
	
	if len(sys.argv)<2:
		sys.exit("  Please use command (python DIANNParamsGen.py d_DIANNParamsOption.params) to re-run")
	params_file = sys.argv[1]
	params = DIANNParamsLib.Get_Params(params_file)
	
	# cur_path
	cur_path = os.getcwd()
	
	# data_path
	data_path = params["data_path"]
	
	# data_path_ln (softlink of data_path)
	# data_path_ln = os.path.join( os.path.split(cur_path)[0],'data')
	# DIANNParamsLib.Check_Path(data_path_ln)
	
	# check data_path (dFolderlist,mzMLlist)
	dFolderlist = DIANNParamsLib.Get_dFolderlist(data_path)
	mzMLlist = DIANNParamsLib.Get_mzMLlist_without_ext(data_path)
	rawlist = DIANNParamsLib.Get_rawlist_without_ext(data_path)
	if len(dFolderlist)==0 and len(mzMLlist)==0 and len(rawlist)==0:
		print('please input the data of .d folders or .mzML files or .raw files in {}'.format(data_path))
		return
	
	# fasta_file
	fasta_file = params["fasta_file"]
	
	# check fasta_file
	fastas = params['fasta_file'].split(";")
	for fasta in fastas:
		if False==os.path.isfile(fasta):
			print('please input the correct .fasta file in {}\ncurrent incorrect .fasta file: {}'.format(os.path.split(fasta)[0], fasta))
			return
	
	# update_date
	update_date = DIANNParamsLib.GetUpdateDate()
	
	
# 1.generate params and sh files #
	if 'Linux'!=platform.system():
		print('please run in Linux')
		return
	
	# cur_path,data_path,fasta_file,dFolderlist,mzMLlist,rawlist,update_date,params
	
	# diann_LibGen_params, diann_Search_params, diann_iqQuan_params
	DIANNParamsLib.generate_diann_LibGen_params(cur_path,data_path,fasta_file,dFolderlist,mzMLlist,rawlist,update_date,params)
	DIANNParamsLib.generate_diann_Search_params(cur_path,data_path,fasta_file,dFolderlist,mzMLlist,rawlist,update_date,params)
	DIANNParamsLib.generate_diann_iqQuan_params(cur_path,data_path,fasta_file,dFolderlist,mzMLlist,rawlist,update_date,params)
	
	# workflow
	DIANNParamsLib.generate_workflow(cur_path,data_path,fasta_file,dFolderlist,mzMLlist,rawlist,update_date,params)
	
	# print
	print('\nDIANN parameters are generated. You can check the parameters and run "bash workflow.sh"\n')
	
# end #
	# end_time = time.time()
	# total_time = end_time-start_time
	# print('\ntotal running time: %.1f sec (%.1f min)' % (total_time,total_time/60))
	return

if __name__ == '__main__':
	main()
