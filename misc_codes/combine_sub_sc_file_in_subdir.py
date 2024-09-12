import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')

def combine_sc_file_in_sub_dir(file_path,ppn=None,file_name=None):
	"""
	This function combines the docking results from sub dir 
	into one giant results file in csv format.

	The result is saved outside the sub dir as docking_result.csv
	Useage : See below after my comment ## This part is my directory specific
	"""
	if file_name == None:
		file_name = "score_file.sc"
	data_frames= []
	if ppn == None:
		ppn = 40
	for i in range(ppn):
		i = i+1
		string_i = str(i)
		df = pd.read_csv(f"{file_path}/{string_i}/{file_name}",sep="\s+",header=1)
		df['Directory']= f"{file_path}/{string_i}/"
		data_frames.append(df)

	combined_score_file = pd.concat(data_frames)
	combined_score_file.to_csv(f"{file_path}/docking_result.csv",index=False)
	i_sc   = combined_score_file["I_sc"]
	f_rmsd = combined_score_file["symmetric_rms"]
	min_rmsd = np.round(min(f_rmsd),2)
	min_isc_idx = i_sc.idxmin()
	min_isc_rmsd= np.round(f_rmsd.iloc[min_isc_idx],2)#np.round(min_isc_idx["rms"],2)

	print(file_path , len(i_sc), min_rmsd, min_isc_rmsd)


## This part is my directory specific ##

af2_protein = pd.read_csv("../homomers/af2_list",header=None)[0]
for p, prot in enumerate(af2_protein):
	file_path_ = f"/fs/project/PAS1146/turzo.1/complex/ideal_complex/2mers/homomers/AF2_DOCK/{prot}"
	#print(file_path_)
	combine_sc_file_in_sub_dir(file_path_,25)

## End of Script ##

