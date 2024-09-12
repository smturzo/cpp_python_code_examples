from pymol import cmd, stored
import pandas as pd

def get_protein_list(file_name):
        with open(file_name) as fileopen:
                lines= fileopen.read().splitlines()
                return lines

def single_chain_rmsd(decoy,native):
	cmd.delete('all')
	cmd.load(str(native), 'native')
	cmd.load(str(decoy),'decoy')
	rmsd_value = float(str(cmd.align('native and n. CA','decoy and n. CA')).split()[3].split(',')[0])
	return [native,decoy,rmsd_value]

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
	print(base_path)
        
	with open(f"{base_path}/{version}_single_chain_subunit_rmsd.csv","w+") as complexfile:
		header_string = "PDBID_CHAINID"
		for af2num in range(num_of_af2_pred):
			header_string = header_string + f",RMSD_R{af2num}"
		complexfile.write(f"{header_string}\n")

		for r, relax_nat in enumerate(relaxed_natives):
			relax_strip_ = relax_nat.split("/")[-1].split(".")[0].split("_")
			decoy_suffix = None
			if relax_strip_[-1] != "A":
				decoy_suffix = "_renamed.pdb"
			else:
				decoy_suffix = ".pdb"
			relax_strip = f"{relax_strip_[1]}_{relax_strip_[-1]}"
			#print(decoy_suffix, relax_strip_[-1])  # Before running RMSD calc. Comment everything below and uncomment the print of this line to make sure it makes sense
			my_rmsd_string =""
			my_rmsd_string = my_rmsd_string + f"{relax_strip}"
			for a in range(num_of_af2_pred):
				decoy_path = f"{base_path}/{version}/{relax_strip}/ranked_{a}{decoy_suffix}"
				if not os.path.exists(decoy_path):
					print(decoy_path,os.path.exists(decoy_path)," This File doesn't exist!, Calculation will fail after this.")
				rmsd_value = single_chain_rmsd(decoy_path,relax_nat)[2]
				string_rmsd_value = str(rmsd_value)
				my_rmsd_string = my_rmsd_string+f",{string_rmsd_value}"

			complexfile.write(f"{my_rmsd_string}\n")

store_af2_rmsd_result_for_complex_decoy("V222","2","homomers","relaxed_single_chain_natives.txt")
store_af2_rmsd_result_for_complex_decoy("V222","2","heteromers","relaxed_single_chain_natives.txt")

store_af2_rmsd_result_for_complex_decoy("V231","2","homomers","relaxed_single_chain_natives.txt")
store_af2_rmsd_result_for_complex_decoy("V231","2","heteromers","relaxed_single_chain_natives.txt")
