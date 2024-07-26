import os, glob
import pandas as pd

# 1
def get_path_of_lowest_energy_structure(directory):
    sc_df = pd.read_csv(f"{directory}/relax_after_genkic_score.sc",sep="\s+",header=1)
    min_score_index = sc_df["total_score"].idxmin()
    lowest_score_struct = sc_df.loc[min_score_index, "description"]
    return f"{directory}/{lowest_score_struct}.pdb"


# 2
def make_slurm_openmm_job(pdbfile,directory,slurm_dir,vac=False):
    if vac == False:
        python_command = f"python /mnt/home/bturzo/ceph/Projects/Holography/04_Scripts/openmm_md_scripts/openmm_driver.py -ipdb {pdbfile} -solvent"
    else:
        python_command = f"python /mnt/home/bturzo/ceph/Projects/Holography/04_Scripts/openmm_md_scripts/openmm_driver.py -ipdb {pdbfile}"
    design_name = pdbfile.split("/")[-1].split(".")[0]
    cluster_slurm_content= f"""#!/bin/bash -l
#SBATCH --job-name={design_name}_MD_{directory.split("/")[-2]}
#SBATCH --nodes=1
#SBATCH -p gpu -C a100
#SBATCH --gpus-per-node=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=140:00:00
#SBATCH -e {directory}/{design_name}_openmm.err
#SBATCH -o {directory}/{design_name}_openmm.log
#SBATCH --mail-user=bturzo@flatironinstitute.org
#SBATCH --mail-type=END,FAIL

module purge
module load modules/2.2
export MODULEPATH=/mnt/home/gkrawezik/modules/rocky8:$MODULEPATH
module load python cuda fftw openmpi openmm/8.1.1
source /mnt/home/gkrawezik/VIRTUAL_ENVIRONMENTS/openmm-8.1.1/bin/activate
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd {directory} 
{python_command}
"""
    with open(f"{slurm_dir}/{design_name}_MD_slurmjob.sh", "w+") as file:
        file.write(cluster_slurm_content)

def print_job_to_run(slurm_dir, design_name,design_flex_type, e_no):
    print(f"cd {slurm_dir} && sbatch {design_name}_MD_slurmjob.sh #{e_no}")

# In case
def clean_out_dir(directory):
    os.system(f"rm {directory}/*.*")

def make_dir(path_to_create_dir,folder_name):
    os.system(f"mkdir -p {path_to_create_dir}/MD/Solution/{folder_name}")


# Base Paths
nstruct_outpath = "/mnt/ceph/users/bturzo/Projects/Holography/08_Protein_Redesign/outputs/"
submit_basepath ="/mnt/ceph/users/bturzo/Projects/Holography/08_Protein_Redesign/submits/"
genkic_design   = pd.read_csv(f"{nstruct_outpath}genkic_designed_proteins.txt",header=0, sep="\s+")
flex_type       = genkic_design["FLEX_TYPE"]
genkic_design_list = genkic_design["GENKIC_DESIGNS"]
ensemble_no        = genkic_design["ENSEMBLE_NUMBER"]

foldr_name = "LowestRosettaEnergy"
#count = 0
append_fix=""
for d, design in enumerate(genkic_design_list):
    genkic_output_dir_= nstruct_outpath+design
    flex_type_ = flex_type[d]
    les = get_path_of_lowest_energy_structure(genkic_output_dir_)
    fix = les.split(".pdb")[0]+"_fix.pdb"
    append_fix +=f"{fix} "
    e_no_ = ensemble_no[d]

    fix_design_name = les.split("/")[-1].split(".pdb")[0]+"_fix"
    genkic_md_slurm_dir_= f"{submit_basepath}/MD/"
    solution_slurm_dir_= f"{genkic_md_slurm_dir_}/Solution"
    vac_slurm_dir_= f"{genkic_md_slurm_dir_}/Vac"
    ## Run Functions
    ##
    solution_md_output_dir = genkic_output_dir_+f"/MD/Solution/{foldr_name}/"
    vac_md_output_dir = genkic_output_dir_+f"/MD/Vac/{foldr_name}/"
    #print(f"rm {vac_md_output_dir}/*.*")
    make_dir(f"{genkic_output_dir_}/",f"{foldr_name}")
    vac = False
    make_slurm_openmm_job(fix, solution_md_output_dir, solution_slurm_dir_, vac)
    print_job_to_run(solution_slurm_dir_, fix_design_name, flex_type_, e_no_)
    
    #vac = True
    #make_slurm_openmm_job(fix, vac_md_output_dir, vac_slurm_dir_, vac)
    #print_job_to_run(vac_slurm_dir_, fix_design_name, flex_type_, e_no_)
#print(append_fix)
