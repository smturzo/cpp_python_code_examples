import pandas as pd
import os
import sys
"""
This script contains function to create a work flow for symmetry docking.
Note in order to do this backbone needs to be done a head of time.
This script assumes that backbone has been moved with Backrub, NMA, Relax and ShearMoves.
The file header of these files are at row 1. So header=1 (for all these 4 files)
Additionally it assumes certain directory structure
For example:
/fs/project/PAS1146/turzo.1/complex/ideal_complex/2mers/homomersrs/XTAL_DOCK/1ECY_AB/
Here for the ideal complex dataset, we are looking at 2mers (only chain A and chain B) 
for the bound crystal structure (XTAL_DOCK) 1ECY_AB. Here 1ECY is the PDB ID and "_AB"
refers to the two chain (A and B), where A and B are identicl thus symmetric dimer or
homodimer.

SCRIPT USAGE on OSC: 
module load python/3.6-conda5.2
python create_symmetry_dock_flags.py > source_this_to_start_hom_2mers_docking_job.sh 
"""


def create_ensemble_list(protein,chain_id,perturbation_file_prefix=None,base_path=None):
	"""
	This function has been designed to read
	/fs/project/PAS1146/turzo.1/complex/ideal_complex/2mers/homomers/XTAL_RELAX_LIST.txt
	Where each line in the file looks like:
	/fs/project/PAS1146/turzo.1/complex/ideal_complex/2mers/homomers/XTAL_DOCK/1ECY_AB/r_1ECY_AB_0001.pdb,C2.symm
	...
	...(EOF)

	Additionally, it assumes that the perturbation score file has the _PerturbationMode_score.sc suffix,
	Where the PerturbationMode currently are nma, relax, backrub and shear.
	for example if you add the prefix full_bound_ and the chain_id=A the function will look for the files
	full_bound_backrub_A_score.sc, full_bound_nma_A_score.sc
	full_bound_relax_A_score.sc, and full_bound_shear_A_score.sc
	"""
	if perturbation_file_prefix == None:
		perturbation_file_prefix = "full_bound_"
	if base_path == None:
		base_path= os.path.dirname(protein)
	ensemble_list_file_name=f"ensemble_list_{chain_id}"
	backrub_sc = pd.read_csv(f"{base_path}/{perturbation_file_prefix}backrub_{chain_id}_score.sc",sep="\s+",header=1)["description"]
	nma_sc   = pd.read_csv(f"{base_path}/{perturbation_file_prefix}nma_{chain_id}_score.sc",sep="\s+",header=1)["description"]
	relax_sc = pd.read_csv(f"{base_path}/{perturbation_file_prefix}relax_{chain_id}_score.sc",sep="\s+",header=1)["description"]
	shear_sc = pd.read_csv(f"{base_path}/{perturbation_file_prefix}shear_{chain_id}_score.sc",sep="\s+",header=1)["description"]
	with open(f"{base_path}/{ensemble_list_file_name}","w+") as ensemble_list_file:
		if len(backrub_sc) == len(nma_sc) == len(relax_sc) == len(shear_sc):
			print("#Equal number of perturbation across all different types of moves used. \n#Will create new file to make the ensemble list")
			for i, bb in enumerate(backrub_sc):
				ensemble_list_file.write(f"{base_path}/{backrub_sc[i]}.pdb\n")
				ensemble_list_file.write(f"{base_path}/{nma_sc[i]}.pdb\n")
				ensemble_list_file.write(f"{base_path}/{relax_sc[i]}.pdb\n")
				ensemble_list_file.write(f"{base_path}/{shear_sc[i]}.pdb\n")
		else:
			print("#MISMATACH in total number of Backbone perturbation. Check your initial perturbatipn score.sc files\n")

def create_dock_dir_and_symmetry_dock_command(protein,chain_A,symm_def,job_name,time="30",ppn=25,decoy_num=10000,rosetta_dock_bin_path=None):
	"""
	This function takes the protein native file with full path
	from which it determines the base path and sets up docking protocols
	It is assuming that the native file exists within the base directory.

	This function assumes you have already done PrePack and 
	the ensemble lists exists accordingly.
	If you need to add or change docking flags, do so in the dock_flags
	variable. 
	"""
	ntsk = str(ppn+2)
	base_path, filename = os.path.split(protein)
	ensemble_list_file_name = f"ensemble_list_{chain_A}"
	perturbed_subunit = pd.read_csv(f"{base_path}/{ensemble_list_file_name}",header=None)[0]
	num_perturbed_subunits= len(perturbed_subunit)
	symm_def_file = f"{base_path}/{symm_def}"
	sub_decoy_num = int((int(decoy_num/ppn))/num_perturbed_subunits)
	if rosetta_dock_bin_path == None:
		rosetta_dock_bin_path = "/fs/project/PAS1146/turzo.1/Rosetta_clone/Rosetta/main/source/bin/SymDock.linuxgccrelease"
	with open(f"{base_path}/dock_commands","w+") as dock_command:
		for i in range(ppn):
			i = i+1
			directory = f"{base_path}/{i}"
			dock_command.write(f"cd {directory}/ && {rosetta_dock_bin_path} @{base_path}/dock_flags\n")
			if os.path.exists(directory):
				print(f"#{i}")
				continue
			else:
				os.mkdir(directory)

	dock_flags = f"""-in:file:l {base_path}/{ensemble_list_file_name}
-in:file:native {protein}
-nstruct {sub_decoy_num}
-symmetry:perturb_rigid_body_dofs 5 60
-symmetry:symmetry_definition {symm_def_file}
-symmetry:initialize_rigid_body_dofs
-symmetry:symmetric_rmsd
-packing:ex1
-packing:ex2aro
-out:prefix dock_
-out:path:all ./
-out:file:silent silent_file.out
-out:file:fullatom
-out:file:scorefile  score_file.sc
-detect_disulf false
-ignore_unrecognized_res

#-evaluation:Irms
#-evaluation:Ca_Irms 
#-evaluation:Fnat 
#-evaluation:Lrmsd 
#-evaluation:Fnonnat 
#-evaluation:contact_map
-evaluation:gdtmm
-evaluation:gdttm

-multiple_processes_writing_to_one_directory

"""
	with open(f"{base_path}/dock_flags","w+") as dock_flg:
		dock_flg.write(dock_flags)
	sub_decoy_num = int((int(decoy_num/ppn))/num_perturbed_subunits)
	sub_decoy_num = int((int(decoy_num/ppn))/num_perturbed_subunits)
	sub_decoy_num = int((int(decoy_num/ppn))/num_perturbed_subunits)
	osc_submit_script = f"""#!/bin/bash
#SBATCH --account=PAS1146
#SBATCH --job-name={job_name}
#SBATCH --time={time}:00:00
#SBATCH --nodes=1 --ntasks={ntsk}
#SBATCH --error=dock_error.txt
#SBATCH --output=dock_log.txt
#SBATCH --mail-user=turzo.1@osu.edu
#SBATCH --mail-type=END,FAIL 

set -vx
module load pcp/pcp
srun parallel-command-processor dock_commands
"""
	with open(f"{base_path}/osc_dock_submit.sh", "w+") as osc_file:
		osc_file.write(osc_submit_script)

	print(f"cd {base_path}/ && sbatch osc_dock_submit.sh")		

# WARNING HARD-CODED DIRECTORY BELOW:
# WARNING WARNING: DO I HAVE YOUR ATTENTION? HARD-CODED DIRECTORY BELOW.
# User can ARG PARSE this. Too Lazy to do this myself.

af2_full_path_list = pd.read_csv("/fs/project/PAS1146/turzo.1/complex/ideal_complex/2mers/homomers/AF2_RELAX_LIST.txt",header=None)
af2_all_relaxed_native = af2_full_path_list[0]
af2_all_symm_def_files = af2_full_path_list[1]

for i, af2 in enumerate(af2_all_relaxed_native):
	relaxed_native_fp = af2_all_relaxed_native[i]
	symm_def_file = af2_all_symm_def_files[i]
	print(f"#{relaxed_native_fp}")
	# Step 1:
	# Use the create ensemble function to create the list of ensemble
	# Use the function twice for chain_A
	create_chainid_A_ensemble = create_ensemble_list(relaxed_native_fp,"A","af2_subunit_v222_")
	# Step 2:
	# Use the create dock flags and osc scripts to create flags and commands to run symmetry docking
	osu_job_jame = relaxed_native_fp.split("/")[-1].split(".")[0]+"_dock_hom_2mers_af2"
	create_dock_dir_and_symmetry_dock_command(relaxed_native_fp,"A",symm_def_file,osu_job_jame)



