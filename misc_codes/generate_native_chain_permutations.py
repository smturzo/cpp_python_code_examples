from create_all_possible_chains import *
import pandas as pd, os, numpy as np, string
from itertools import permutations
import threading
#../homomers/2mers/native_list
def create_dir_for_list_of_natives(native):	
	protein_dir = os.path.abspath(os.path.dirname(native))
	#print(protein_dir)
	os.system("mkdir -p "+protein_dir+"/CNAT")
	base_dir = os.path.dirname(os.path.dirname(protein_dir))
	protein_name = "/r_"+protein_dir.split("/")[-1]+"_0001.pdb"
	#print(base_dir+protein_name)
	output_nat_dir = protein_dir+"/CNAT/"
	print("ls "+output_nat_dir+"*_*.pdb > "+output_nat_dir+"gen_list")
	generate_all_chain_combination_for_pdb(base_dir+protein_name,output_nat_dir)

native_list= pd.read_csv("../heteromers/v231_complex_AB_ranked_0_list.txt",header=None)[0]
for i, nat in enumerate(native_list):
	my_thread = threading.Thread(target=create_dir_for_list_of_natives,args=(nat,))
	my_thread.start()

