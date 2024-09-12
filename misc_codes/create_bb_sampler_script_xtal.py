import numpy as np
import pandas as pd
def full_backbone_move(protein,complex_type,mer,chain):
    rosetta_bin_path ="/fs/project/PAS1146/turzo.1/Rosetta_clone/Rosetta/main/source/bin/"
    rosetta_db_path  ="-database /fs/project/PAS1146/turzo.1/Rosetta_clone/Rosetta/main/database"
    protein_path=f"/fs/project/PAS1146/turzo.1/complex/ideal_complex/{mer}/{complex_type}/XTAL_DOCK/{protein}/"
    native_path=f"{protein_path}r_{protein}_0001_{chain}.pdb"

    """ Create the ensemble list for backbone perturbation """
    print(f"cd {protein_path} && sbatch full_bound_crystal_backbone_move_job_{chain}.sh")
    with open(f"{protein_path}full_bound_crystal_backbone_move_job_{chain}.sh","w+") as full_crystal_osc_file:
        bound_commands_to_write= f'#!/bin/bash\n#SBATCH --account=PAS1146\n#SBATCH --job-name={protein}_{chain}_{complex_type}_{mer}\n#SBATCH --time=60:00:00\n#SBATCH --nodes=1 --ntasks=20\n#SBATCH --error=full_perturb_error.txt\n#SBATCH --output=full_perturb_log.txt\n#SBATCH --mail-user=turzo.1@osu.edu\n#SBATCH --mail-type=END,FAIL\nset -vx\nmodule load pcp/pcp\nsrun parallel-command-processor full_bound_bb_starter_{chain}\n'
        full_crystal_osc_file.write(bound_commands_to_write)
    with open(f"{protein_path}full_bound_bb_starter_{chain}",'w+') as bound_bb_file:
        bound_bb_file.write(f"cd {protein_path} && {rosetta_bin_path}backrub.default.linuxgccrelease {rosetta_db_path} -in:file:s {native_path} -nstruct 10  -backrub:mc_kt 0.6 -backrub:ntrials 200 -out:prefix full_bound_backrub_{chain}_ \n")
        bound_bb_file.write(f"cd {protein_path} && {rosetta_bin_path}relax.linuxgccrelease {rosetta_db_path} -in:file:s {native_path} -relax:quick -nstruct 10 -out:prefix full_bound_relax_{chain}_ \n")
        bound_bb_file.write(f"cd {protein_path} && {rosetta_bin_path}rosetta_scripts.linuxgccrelease {rosetta_db_path} -in:file:s {native_path}  -nstruct 10 -parser:protocol /fs/project/PAS1146/turzo.1/complex/scripts/bb_sampler/full_nma.xml -out:prefix full_bound_nma_{chain}_\n")
        bound_bb_file.write(f"cd {protein_path} && {rosetta_bin_path}rosetta_scripts.linuxgccrelease {rosetta_db_path} -in:file:s {native_path} -nstruct 10 -parser:protocol /fs/project/PAS1146/turzo.1/complex/scripts/bb_sampler/full_shear.xml -ex1 -ex2aro -out:prefix full_bound_shear_{chain}_\n")



xtal_het_data = pd.read_csv("/fs/project/PAS1146/turzo.1/complex/ideal_complex/2mers/heteromers/xtal_list",header=None)[0]
xtal_hom_data = pd.read_csv("/fs/project/PAS1146/turzo.1/complex/ideal_complex/2mers/homomers/xtal_list",header=None)[0]

#for x, xtal in enumerate(xtal_het_data):
#    #print(xtal)
#    full_backbone_move(xtal,"heteromers","2mers","A")
#    full_backbone_move(xtal,"heteromers","2mers","B")

for y, xtall in enumerate(xtal_hom_data):
    full_backbone_move(xtall,"homomers","2mers","A")





