import numpy as np
import os
import pandas as pd

def equil_slurm_config(pdb_id, output_dir, init="step3_input"):
    slurm_equil = f"""#!/bin/bash
#SBATCH --job-name=Equil_{pdb_id}
#SBATCH --nodes=1
#SBATCH -p gpu -C a100
#SBATCH --gpus-per-node=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=100:00:00
#SBATCH -e {output_dir}/err.equil.{pdb_id}
#SBATCH -o {output_dir}/log.equil.{pdb_id}
#SBATCH --mail-user=bturzo@flatironinstitute.org
#SBATCH --mail-type=END,FAIL

module purge
module load modules/2.2
#export MODULEPATH=/mnt/home/gkrawezik/modules/rocky8:$MODULEPATH
#module load python cuda fftw openmpi openmm/8.1.1
#source /mnt/home/gkrawezik/VIRTUAL_ENVIRONMENTS/openmm-8.1.1/bin/activate
#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module use /mnt/ceph/users/protfold/bmd_modules/modules/
module load open_mm/8.1.0.beta


#mkdir -p {output_dir}
cd {output_dir}

# Equilibration
python /mnt/home/bturzo/ceph/Projects/KRAS_ENSEMBLE_METRICS/scripts/openmm_tremd.py -run_nvt_sim -topparstr ./toppar.str -path_to_psf ./{init}.psf -path_to_crdfile ./{init}.crd 
"""
    with open(f"/mnt/home/bturzo/ceph/Projects/KRAS_ENSEMBLE_METRICS/submits/Submit_Equil_{pdb_id}.sh", "w") as sbmt_equil:
        sbmt_equil.write(slurm_equil)

def trex_slurm_config(pdb_id, output_dir):
    slurm_prod = f"""#!/bin/bash
#SBATCH --job-name=Prod_{pdb_id}
#SBATCH --nodes=1
#SBATCH -p gpu -C a100
#SBATCH --gpus-per-node=1

#SBATCH --cpus-per-task=20
#SBATCH --time=168:00:00
#SBATCH -e {output_dir}/err.trex.{pdb_id}
#SBATCH -o {output_dir}/log.trex.{pdb_id}
#SBATCH --mail-user=bturzo@flatironinstitute.org
#SBATCH --mail-type=END,FAIL

module purge
module load modules/2.2
export MODULEPATH=/mnt/home/gkrawezik/modules/rocky8:$MODULEPATH
module load python cuda fftw openmpi openmm/8.1.1
source /mnt/home/gkrawezik/VIRTUAL_ENVIRONMENTS/openmm-8.1.1/bin/activate
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#module use /mnt/ceph/users/protfold/bmd_modules/modules/
#module load open_mm/8.1.0.beta


#mkdir -p {output_dir}
cd {output_dir}

# Temperature Replica Exchange Production Run 
python /mnt/home/bturzo/ceph/Projects/KRAS_ENSEMBLE_METRICS/scripts/openmm_tremd.py -restart -irstno 1 -run_trex_sim 
"""
    with open(f"/mnt/home/bturzo/ceph/Projects/KRAS_ENSEMBLE_METRICS/submits/trex/restart_1/Submit_Trex_{pdb_id}_Restart_1.sh", "w") as sbmt_prod:
        sbmt_prod.write(slurm_prod)

toppar_files = pd.read_csv("/mnt/home/bturzo/ceph/Projects/KRAS_ENSEMBLE_METRICS/inputs/toppar_list.txt",header=None,sep="\s+")[0]
psf_files    = pd.read_csv("/mnt/home/bturzo/ceph/Projects/KRAS_ENSEMBLE_METRICS/inputs/psf_files.txt",header=None,sep="\s+")[0]
crd_files    = pd.read_csv("/mnt/home/bturzo/ceph/Projects/KRAS_ENSEMBLE_METRICS/inputs/crd_files.txt",header=None,sep="\s+")[0]

for f, fls in enumerate(toppar_files):
    pdbid = crd_files[f].split('/')[-3]
    outputdir = os.path.dirname(psf_files[f])
    print(pdbid, outputdir)
    #equil_slurm_config(pdbid, outputdir)
    trex_slurm_config(pdbid, outputdir)
