import pandas as pd, os, numpy as np, string
from itertools import permutations
#import os.path
#from data import *
def generate_all_chain_combination_for_pdb(file_name,path_to_save):
	alphabet_dict = dict(enumerate(string.ascii_uppercase)) # a dictinary that converts number to uppercase alphabets, think this is pretty cool
	chain_list_from_pdb =[]
	pdb_items_before_chainid=[]
	pdb_items_after_chainid =[]
	pdb_items_with_chainid =[]
	with open(str(file_name),'r') as pdbfile:
		for i, line in enumerate(pdbfile):
			#print(line)

			if line.startswith("ATOM") and line.split()[0]!='TER':	
				pdb_items_before_chainid.append(line[0:21])
				pdb_items_after_chainid.append(line[22:])
				pdb_items_with_chainid.append(line[21])
				if line[21] not in chain_list_from_pdb:
					chain_list_from_pdb.append(line[21])

	#print(len(chain_list_from_pdb),len(pdb_items_after_chainid),len(pdb_items_before_chainid)) #Sanity check
	modified_chain_list =[alphabet_dict.get(x) for x in range(len(chain_list_from_pdb))] #This converts numbers to alphabet based on the number of chanis from og PDB
	#print(modified_chain_list)
	#print(chain_list_from_pdb)
	all_permutation_list = list(permutations(modified_chain_list))
	#print(all_permutation_list)
	perm_made_pdb_list =[]
	for p, perm in enumerate(all_permutation_list):
		permuted_list=list(perm)
			
		pdb_to_modified_chain_dict ={} # Calling a dictionary here so that it can map og pdb chain list to permuted chain list. This needs to be done for every permutation
		for i, chain in enumerate(chain_list_from_pdb):
			pdb_to_modified_chain_dict[chain_list_from_pdb[i]]=permuted_list[i]
		native_prefix =''.join(permuted_list)
		#print(pdb_to_modified_chain_dict)
		perm_made_pdb_list.append(str(path_to_save)+str(file_name).split('/')[-1].split(".pdb")[0]+'_'+str(native_prefix)+'.pdb')
		with open(str(path_to_save)+str(file_name).split('/')[-1].split(".pdb")[0]+'_'+str(native_prefix)+'.pdb','w+') as modified_pdb:
			for i, line in enumerate(pdb_items_before_chainid): 
				modified_pdb.write(str(pdb_items_before_chainid[i])+str(pdb_to_modified_chain_dict.get(pdb_items_with_chainid[i]))+str(pdb_items_after_chainid[i]))
			modified_pdb.write('TER')	

	return perm_made_pdb_list

def get_native_list(user_path,mer,protein):
        gen_path = user_path+mer+"/"+protein+"CNAT/"
        generated_native_pdb_list = []
        with open(gen_path+'gen_list','r') as read_gen_nat_file:
                for i, line in enumerate(read_gen_nat_file):
                        generated_native_pdb_list.append(line[:-1])

        return generated_native_pdb_list



#file_name='3INS_AB.pdb'
#path_to_save='./'
#generate_all_chain_combination_for_pdb(file_name,path_to_save)

#protein_list = protein_list_data()

#two_mers  = protein_list[0]
#four_mers = protein_list[1]
#five_mers = protein_list[2]
#six_mers  = protein_list[3]

#for key in six_mers:
#	file_name="/fs/project/PAS1146/turzo.1/complex/"+six_mers.get(key)+"/"+key+"/AF/HM/"+key[0:4]+".pdb"
#	print(file_name)
#	path_to_save="/fs/project/PAS1146/turzo.1/complex/"+six_mers.get(key)+"/"+key+"/AF/CNAT/"
#	generate_all_chain_combination_for_pdb(file_name,path_to_save)

