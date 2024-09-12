from pymol import cmd, stored
from itertools import permutations
import os

# Chain list obtained from both native and decoy
# within the functions
# This is to handle cases such as if native is A,C, and D
# But decoy has chains A, B, and C. Then this function will
# change the decoy to have chains as the native in that order.

def min_rmsd_perm_any_chain(decoy_path,native_path):
	cmd.delete('all')                                          # Clean Session
	cmd.load(native_path,'native')                              # Load Native
	native_chains = list(cmd.get_chains('native'))             # Get chain list of natives, this is needed for permutations
	num_chains = len(native_chains)                            # Get number of chains in from chain list
	permlist = list(permutations(native_chains, num_chains))   # Create a permutation list (from native but doesn't matter even if from decoy)

	cmd.delete('all')                                          # Clean session again
	cmd.load(decoy_path, 'current')                            # Load Decoy
	decoy_chains = list(cmd.get_chains('current'))             # Get the Decoy chain list, this is needed for altering decoy chain id to temp ids
	cmd.delete('all')                                          # Clean session
	
	min_rmsd = float('inf')  # Set min RMSD to a riddiculously high value as the initial value (since we will minimize this)
	comb_best = ''                                             # Empty string to store the best combo of chain list
	for comb in permlist:                                      # for each combo in the permutations
		cmd.delete('all')                                  # Clean session
		cmd.load(native_path,'native')                     # Load Native
		cmd.load(decoy_path, 'current')                    # Load Decoy
		#decoy_chains = cmd.get_chains('current')        
		for i in range(num_chains):                        # Alter chain of decoy (from decoy chain list) to a temp chain ID
			cmd.alter(f'current and chain {decoy_chains[i]}', f'chain = "tem{i+1}"')
		
		for i in range(num_chains):                        # Alter chain of temp chain ID to the combo
			cmd.alter(f'current and chain tem{i+1}', f'chain = "{comb[i]}"')

		cmd.sort('current')                                # Sort the chains after altering with combo
		align_rmsd = float(cmd.align('current and n. CA','native and n. CA')[3]) # Find RMSD 
		if align_rmsd < min_rmsd:                          # If align rmsd is less than min_rmsd, then update min_rmsd
			min_rmsd = align_rmsd                      # By setting min_rmsd  = align rmsd
			#This list contains the what each chain should be relabeled to in order to get the minimum RMSD
			comb_best = comb                           # Store the best combo
	comb_best_list = list(comb_best)                           # Get the best combo as a list
	cmd.delete('all')                                          # Clean Session again
	cmd.load(decoy_path,'current')                             # Load the Decoy
	for i in range(num_chains):                                # Alter decoy chain to temp ID
		cmd.alter(f'current and chain {decoy_chains[i]}', f'chain = "tem{i+1}"')
	for i in range(num_chains):                                # Alter temp chain ID to best combo chain ID
		cmd.alter(f'current and chain tem{i+1}', f'chain = "{comb_best_list[i]}"')
	cmd.sort('current')                                        # Sort and this should now correspond to the native
	base_path, filename = os.path.split(decoy_path)
	new_filename = filename.split(".")[0]
	cmd.save(f'{base_path}/{new_filename}_renamed.pdb','current')# Save decoy.
	return [native_path,decoy_path,min_rmsd]


def store_af2_rmsd_result_for_complex_decoy(version,nsubunits,complex_type,relax_native_filename,num_of_af2_pred=5,user_defined_relaxed_native_filepath=None):
	"""
	This function takes version number in the variable version as either V222 or V231
	This function takes number of total subunits of complexes as a string number "2" not 2
	This function defines complex_type as either heteromers or homomers
	This function takes relax_native_filename as a string. 
	Note that the base path for relax_native_filename is hardcoded to /fs/project/PAS1146/turzo.1/complex/ideal_complex/
	User may chose to add their path to native relaxed file with full path (REQUIRED FOR OTHER USERS NOT ME)
	num_of_af2_pred is the number of AF2 prediction to calculate RMSD for.
	"""
	import pandas as pd # Some pymol versions does not support pandas
	if user_defined_relaxed_native_filepath == None:
		relax_native_filenames = f"/fs/project/PAS1146/turzo.1/complex/ideal_complex/{nsubunits}mers/{complex_type}/{relax_native_filename}"
	else:
		relax_native_filenames = f"{user_defined_relaxed_native_filepath}/{relax_native_filename}"

	base_path, filename = os.path.split(relax_native_filenames)
	protein_df = pd.read_csv(relax_native_filenames,header=None,sep=",")
	relaxed_natives = protein_df[0]
	with open(f"{base_path}/{version}_full_complex_rmsd.csv","w+") as complexfile:
		header_string = "PDBID_CHAINID"
		for af2num in range(num_of_af2_pred):
			header_string = header_string + f",RMSD_R{af2num}"
		complexfile.write(f"{header_string}\n")

		for r, relax_nat in enumerate(relaxed_natives):
			relax_strip = relax_nat.split("/")[-1].replace(r'r_', '')
			relax_strip = relax_strip.replace(r'_0001.pdb', '')
			my_rmsd_string =""
			print(relax_strip)
			my_rmsd_string = my_rmsd_string + f"{relax_strip}"
			for a in range(num_of_af2_pred):
				decoy_path = f"{base_path}/{version}/{relax_strip}/ranked_{a}.pdb"
				rmsd_value = min_rmsd_perm_any_chain(decoy_path,relax_nat)[2]
				string_rmsd_value = str(rmsd_value)
				my_rmsd_string = my_rmsd_string+f",{string_rmsd_value}"
			complexfile.write(f"{my_rmsd_string}\n")


store_af2_rmsd_result_for_complex_decoy("V222","2","homomers","relaxed_AB_natives.txt")
store_af2_rmsd_result_for_complex_decoy("V222","2","heteromers","relaxed_AB_natives.txt")

store_af2_rmsd_result_for_complex_decoy("V231","2","homomers","relaxed_AB_natives.txt")
store_af2_rmsd_result_for_complex_decoy("V231","2","heteromers","relaxed_AB_natives.txt")






