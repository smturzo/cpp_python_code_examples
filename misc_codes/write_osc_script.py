import pandas as pd, os

def generate_ascend_af2_script(job_name, fasta_paths, time, module_version, model_preset, max_template_date, output_dir):
     with open(f"{output_dir}AF2_SUBMIT_SCRIPTS/{job_name}_{module_version}_AF2_OSC_script.sh", "w") as f:
          f.write(f"#!/bin/bash\n")
          f.write(f"#SBATCH --account=PAS1146\n")
          f.write(f"#SBATCH --job-name={job_name}\n")
          f.write(f"#SBATCH --time={time}\n")
          f.write(f"#SBATCH --nodes=1\n")
          f.write(f"#SBATCH --ntasks=48\n")
          f.write(f"#SBATCH --gpus-per-node=4\n")
          f.write(f"#SBATCH --gpu_cmode=shared\n")
          f.write(f"#SBATCH --error={job_name}_error.txt\n")
          f.write(f"#SBATCH --output={job_name}_log.txt\n")
          f.write(f"#SBATCH --mail-user=turzo.1@osu.edu\n")
          f.write(f"#SBATCH --mail-type=END,FAIL\n\n")
          f.write(f"set -vx\n")
          f.write(f"module reset\n")
          f.write(f"module load alphafold/{module_version}\n\n")
          if module_version == '2.2.2':
               f.write(f"export TF_FORCE_UNIFIED_MEMORY=1\n")
               f.write(f"run_alphafold.sh ")
               f.write(f"--output_dir={output_dir}V222/ ")
          else:
               f.write(f"#export TF_FORCE_UNIFIED_MEMORY=1\n")
               f.write(f"/fs/project/PAS1146/turzo.1/AF231/alphafold/run_alphafold.sh ")
               f.write(f"--output_dir={output_dir}V231/  ")
          f.write(f"--fasta_paths={fasta_paths}{job_name}.fasta ")
          f.write(f"--db_preset=full_dbs ")
          f.write(f"--model_preset={model_preset} ")
          f.write(f"--max_template_date={max_template_date} ")
          if model_preset == "multimer":
               f.write(f"--num_multimer_predictions_per_model=1 ")
          f.write(f"--use_gpu_relax ")


output_dir="/fs/project/PAS1146/turzo.1/complex/ideal_complex/2mers/heteromers/"
time="10:00:00"
module_version ="2.2.2"
model_preset= "multimer"
max_template_date ="1900-01-01"
fasta_paths=output_dir

homomer_df = pd.read_csv("../heteromers/list_2mers_fasta.txt", header=None)[0].apply(os.path.basename)

for homomer in homomer_df:
    job_name = homomer.split(".")[0]
    generate_ascend_af2_script(job_name, fasta_paths, time, module_version, model_preset, max_template_date, output_dir)
    print("sbatch "+job_name+"_"+module_version+"_AF2_OSC_script.sh")


output_dir="/fs/project/PAS1146/turzo.1/complex/ideal_complex/2mers/heteromers/"
time="06:00:00"
module_version ="2.2.2"
model_preset= "monomer"
max_template_date ="1900-01-01"
fasta_paths=output_dir

homomer_df = pd.read_csv("../heteromers/list_2mers_subunit_fasta.txt", header=None)[0].apply(os.path.basename)

for homomer in homomer_df:
    job_name = homomer.split(".")[0]
    generate_ascend_af2_script(job_name, fasta_paths, time, module_version, model_preset, max_template_date, output_dir)
    print("sbatch "+job_name+"_"+module_version+"_AF2_OSC_script.sh")






