from Bio.PDB import PDBList, PDBIO, PDBParser
import os
import pandas as pd
io = PDBIO()
parser = PDBParser()

def rename_single_chain(af_struct_file):
    io = PDBIO()
    parser = PDBParser()
    structure = parser.get_structure('decoy', str(af_struct_file))
    renames = {"A": "B"}
    for model in structure:
        for chain in model:
            old_name = chain.get_id()
            new_name = renames.get(old_name)
            if new_name:
                print(f"renaming chain {old_name} to {new_name} of file {af_struct_file}")
                chain.id = new_name
            else:
                print(f"keeping chain name {old_name}")
    base_dir_path = os.path.dirname(af_struct_file)
    out_file_prepend = af_struct_file.split("/")[-1].split(".")[0]
    io.set_structure(structure)
    io.save('{}/{}_renamed.pdb'.format(base_dir_path,out_file_prepend))
#
#file_df = pd.read_csv('/fs/project/PAS1146/turzo.1/complex/ideal_complex/2mers/heteromers/all_v222_subunit_b_list.txt',header=None)[0]
file_df = pd.read_csv('/fs/project/PAS1146/turzo.1/complex/ideal_complex/2mers/heteromers/all_v231_subunit_b_list.txt',header=None)[0]

for f_df in file_df:
    #print(f_df)
    rename_single_chain(f_df)

