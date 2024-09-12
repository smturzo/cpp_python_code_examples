import os
import numpy as np
import pandas as pd
def af2_subunit_backbone_move(protein,complex_type,mer,chain,decoy_index,af2_version_number):
    rosetta_bin_path ="/fs/project/PAS1146/turzo.1/Rosetta_clone/Rosetta/main/source/bin/"
    rosetta_db_path  ="-database /fs/project/PAS1146/turzo.1/Rosetta_clone/Rosetta/main/database"
    protein_path=f"/fs/project/PAS1146/turzo.1/complex/ideal_complex/{mer}/{complex_type}/AF2_DOCK/{protein}/"
    decoy_name = f"ranked_{decoy_index}.pdb"
    subunit_name= protein.split("_")[0]+f"_{chain}"
    if chain != "A":
        decoy_name = f"ranked_{decoy_index}_renamed.pdb"
    af2_subunit_path=f"/fs/project/PAS1146/turzo.1/complex/ideal_complex/{mer}/{complex_type}/V{af2_version_number}/{subunit_name}"
    native_path=f"{protein_path}r_{protein}_0001_{chain}.pdb"
    decoy_subunit_path=f"{af2_subunit_path}/{decoy_name}"

    """ Create the ensemble list for backbone perturbation """
    print(f"cd {protein_path} && sbatch af2_subunit_backbone_move_job_{chain}_v{af2_version_number}.sh")
    with open(f"{protein_path}af2_subunit_backbone_move_job_{chain}_v{af2_version_number}.sh","w+") as af2_subunit_osc_file:
        bound_commands_to_write= f'#!/bin/bash\n#SBATCH --account=PAS1146\n#SBATCH --job-name=af2_{subunit_name}_{complex_type}_{mer}_v{af2_version_number}\n#SBATCH --time=10:00:00\n#SBATCH --nodes=1 --ntasks=20\n#SBATCH --error=af2_subunit_perturb_error.txt\n#SBATCH --output=af2_subunit_perturb_log.txt\n#SBATCH --mail-user=turzo.1@osu.edu\n#SBATCH --mail-type=END,FAIL\nset -vx\nmodule load pcp/pcp\nsrun parallel-command-processor af2_subunit_bb_starter_{chain}\n'
        af2_subunit_osc_file.write(bound_commands_to_write)
    with open(f"{protein_path}af2_subunit_bb_starter_{chain}",'w+') as bound_bb_file:
        bound_bb_file.write(f"cd {protein_path} && {rosetta_bin_path}backrub.default.linuxgccrelease {rosetta_db_path} -in:file:s {decoy_subunit_path} -nstruct 10  -backrub:mc_kt 0.6 -backrub:ntrials 200 -out:prefix af2_subunit_v{af2_version_number}_backrub_{chain}_ \n")
        bound_bb_file.write(f"cd {protein_path} && {rosetta_bin_path}relax.linuxgccrelease {rosetta_db_path} -in:file:s {decoy_subunit_path} -relax:quick -nstruct 10 -out:prefix af2_subunit_v{af2_version_number}_relax_{chain}_ \n")
        bound_bb_file.write(f"cd {protein_path} && {rosetta_bin_path}rosetta_scripts.linuxgccrelease {rosetta_db_path} -in:file:s {decoy_subunit_path}  -nstruct 10 -parser:protocol /fs/project/PAS1146/turzo.1/complex/scripts/bb_sampler/full_nma.xml -out:prefix af2_subunit_v{af2_version_number}_nma_{chain}_\n")
        bound_bb_file.write(f"cd {protein_path} && {rosetta_bin_path}rosetta_scripts.linuxgccrelease {rosetta_db_path} -in:file:s {decoy_subunit_path} -nstruct 10 -parser:protocol /fs/project/PAS1146/turzo.1/complex/scripts/bb_sampler/full_shear.xml -ex1 -ex2aro -out:prefix af2_subunit_v{af2_version_number}_shear_{chain}_\n")



af2_het_data = pd.read_csv("/fs/project/PAS1146/turzo.1/complex/ideal_complex/2mers/heteromers/af2_list",header=None)[0]
af2_hom_data = pd.read_csv("/fs/project/PAS1146/turzo.1/complex/ideal_complex/2mers/homomers/af2_list",header=None)[0]

for x, af2sub in enumerate(af2_hom_data):
    af2_subunit_backbone_move(af2sub,"homomers","2mers","A","0","222")


for y, af2sub in enumerate(af2_het_data):
    af2_subunit_backbone_move(af2sub,"heteromers","2mers","A","0","222")
    af2_subunit_backbone_move(af2sub,"heteromers","2mers","B","0","222")





