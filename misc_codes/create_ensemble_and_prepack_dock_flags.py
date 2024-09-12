import pandas as pd
import os
import sys
"""
This script contains function to create a work flow for ensemble docking.
Note in order to do this backbone needs to be done a head of time.
This script assumes that backbone has been moved with Backrub, NMA, Relax and ShearMoves.
The file header of these files are at row 1. So header=1 (for all these 4 files)
Additionally it assumes certain directory structure
For example:
/fs/project/PAS1146/turzo.1/complex/ideal_complex/2mers/heteromers/XTAL_DOCK/6I82_AB/
Here for the ideal complex dataset, we are looking at 2mers (only chain A and chain B) 
for the bound crystal structure (XTAL_DOCK) 6I82_AB. Here 6I82 is the PDB ID and "_AB"
refers to the two chain (A and B).

SCRIPT USAGE on OSC: 
module load python/3.6-conda5.2
python create_ensemble_and_prepack_dock_flags.py > source_this_to_start_het_2mers_docking_job.sh 
"""


def create_ensemble_list(protein,chain_id,perturbation_file_prefix=None,base_path=None):
	"""
	This function has been designed to read
	/fs/project/PAS1146/turzo.1/complex/ideal_complex/2mers/heteromers/XTAL_RELAX_LIST.txt
	Where each line in the file looks like:
	/fs/project/PAS1146/turzo.1/complex/ideal_complex/2mers/heteromers/XTAL_DOCK/1EX4_AB/r_1EX4_AB_0001.pdb
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

def create_and_run_prepack_flags(protein,native,chain_A,chain_B,ppk_prefix="ppk_ensemble_",rosetta_ppk_bin_path=None):
	"""
	This function creates the flags and commands to run prepack for ensemble docking
	The function takes the native protein with full path and obtains the base path. 
	The function also takes chain ids of the docking partners "A" and "B" for a dimeric complex.
	Note this function assumes that Rosetta Database
	can be found in path
	"""
	
	if rosetta_ppk_bin_path == None:
		rosetta_ppk_bin_path = "/fs/project/PAS1146/turzo.1/Rosetta_clone/Rosetta/main/source/bin/docking_prepack_protocol.linuxgccrelease"
	base_path=os.path.dirname(protein)
	ppk_job_name = base_path.split("/")[-1]
	flag_file_name = "ppk_flags"
	ppk_flags=f"""\
-in:file:s {protein}
-in:file:native {native}

-nstruct 1
-partners {chain_A}_{chain_B}
-ensemble1 {base_path}/ensemble_list_{chain_A}
-ensemble2 {base_path}/ensemble_list_{chain_B}

-ex1
-ex2aro
-out:prefix {ppk_prefix}"""
	#print(ppk_flags)
	with open(f"{base_path}/{flag_file_name}","w+") as ppk_f:
		ppk_f.write(ppk_flags)
	with open(f"{base_path}/run_ppk.sh","w+") as run_ppk:

		run_ppk.write(f"#!/bin/bash\n#SBATCH --account=PAS1146\n#SBATCH --job-name={ppk_prefix}_{ppk_job_name}\n")
		run_ppk.write(f"##SBATCH --time=09:00:00\n")
		run_ppk.write(f"#SBATCH --nodes=1 --ntasks=42\n")
		run_ppk.write(f"#SBATCH --error=ppk_error.txt\n")
		run_ppk.write(f"#SBATCH --output=ppk_log.txt\n")
		run_ppk.write(f"#SBATCH --mail-user=turzo.1@osu.edu\n")
		run_ppk.write("#SBATCH --mail-type=END,FAIL\nset -vx\n")
		run_ppk.write(f"cd {base_path} && {rosetta_ppk_bin_path} @{flag_file_name}")
	print(f"#cd {base_path} && sbatch run_ppk.sh")

def create_dock_dir_and_ensemble_dock_command(protein,chain_A,chain_B,job_name,time="20",ppk_prefix="ppk_ensemble_",ppn=40,decoy_num=10000,rosetta_dock_bin_path=None):
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
	sub_decoy_num = int(decoy_num/ppn)
	base_path, filename = os.path.split(protein)
	temp_filename = filename.split(".")[0]
	new_filename = f"{ppk_prefix}{temp_filename}_0001.pdb"
	ppk_prot = os.path.join(base_path, new_filename)
	if rosetta_dock_bin_path == None:
		rosetta_dock_bin_path = "/fs/project/PAS1146/turzo.1/Rosetta_clone/Rosetta/main/source/bin/docking_protocol.linuxgccrelease"
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

	dock_flags = f"""-in:file:s {ppk_prot} 
-in:file:native {protein}
-nstruct {sub_decoy_num}
-partners {chain_A}_{chain_B}
-randomize2

-ensemble1 {base_path}/ensemble_list_{chain_A}.ensemble
-ensemble2 {base_path}/ensemble_list_{chain_B}.ensemble
-ex1
-ex2aro
-out:prefix dock_
-out:path:all ./
-out:file:silent silent_file.out
-out:file:fullatom
-out:file:scorefile  score_file.sc
-detect_disulf false
-ignore_unrecognized_res

-evaluation:Irms
-evaluation:Ca_Irms 
-evaluation:Fnat 
-evaluation:Lrmsd 
-evaluation:Fnonnat 
-evaluation:contact_map
-evaluation:gdtmm
-evaluation:gdttm

-multiple_processes_writing_to_one_directory
"""
	with open(f"{base_path}/dock_flags","w+") as dock_flg:
		dock_flg.write(dock_flags)
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

af2_full_path_list = pd.read_csv("/fs/project/PAS1146/turzo.1/complex/ideal_complex/2mers/heteromers/AF2_RELAX_LIST.txt",header=None)[0]

for i, af2 in enumerate(af2_full_path_list):
	print(f"#{af2}")
	# Step 1:
	# Use the create ensemble function to create the list of ensemble
	# Use the function twice for chain_A and chain_B
	create_chainid_A_ensemble = create_ensemble_list(af2,"A","af2_subunit_v222_")
	create_chainid_B_ensemble = create_ensemble_list(af2,"B","af2_subunit_v222_")
	# Step 2:
	# Use the create pre_pak flags function to generate the flags and run commands for
	# running prepack for every protein in the loop. This function can take both chains and you
	# only need to run it once
	create_and_run_prepack_flags(af2,af2,"A","B")
	print("\n\n#RUN PPK BEFORE RUNNING DOCKING!!\n\n")
	# Step 3:
	# Use the create dock flags and osc scripts to create flags and commands to run ensemble docking
	# This function can take both chain_A and chain_B at the same time and only needs to be run once.
	create_dock_dir_and_ensemble_dock_command(af2,"A","B",af2.split("/")[-1].split(".")[0]+"_dock_het_2mers_af2")
