from Bio.PDB import *
import pandas as pd

def get_avg_plddt(protein):
	parser = PDBParser(QUIET = True)
	data = parser.get_structure("PROT", str(protein))
	protein_name= protein.split("/")[-2]
	#print(data,type(data))
	chains = list(data.get_chains())  
	#outf =open(str(outfile),"w+")
	#outf.write("CHAIN_ID,RES_ID,RES_NAME,pLDDT_by_cutoff,\n")
	total_ca_plddt = 0
	total_number_of_residues = 0
	for c, chain in enumerate(chains):
		chain_id = str(chain.get_id())
		residues = list(chain.get_residues())
		total_number_of_residues += len(residues)
		#print("Total number of residues found: "+str( len(residues) )+" in chain "+chain_id+" for protein "+protein)
		for i, res in enumerate(residues):
			atoms = list(res.get_atoms())
			atoms.sort()
			#print(atoms[1])
			total_ca_plddt += atoms[1].get_bfactor()
	avg_plddt = total_ca_plddt / total_number_of_residues
	print(protein_name+" "+str(avg_plddt))

ideal_ds_A = pd.read_csv("../heteromers/v222_subunit_A_ranked_0_list.txt",header=None)[0]
ideal_ds_B = pd.read_csv("../heteromers/v222_subunit_B_ranked_0_list.txt",header=None)[0]
for p, _ in enumerate(ideal_ds_A):
	get_avg_plddt(ideal_ds_A[p])
	get_avg_plddt(ideal_ds_B[p])
