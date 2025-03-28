import sys
import argparse
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
from math import *
import os

def get_water_count(pdb):
    """
    Count the total number of water res in a given PDB structure.

    This function iterates through the residues of the provided PDB structure
    and counts the number of residues that belong to water molecules. It recognizes
    water residues by the names 'HOH' and 'WAT'.

    Parameters:
    ----------
    pdb : PDBFile
        A PDBFile object containing the molecular structure with topology and atom information.

    Returns:
    -------
    water_atom_count : int
        The total number of water residues in the PDB structure.

    Notes:
    ------
    - The function assumes that the input PDB file has already been parsed into a PDBFile object.
    - Residues with names 'HOH' or 'WAT' are considered water molecules.

    Examples:
    --------
    # pdb = PDBFile('water_structure.pdb')
    # water_count = get_water_count(pdb)
    # print(f"Number of water atoms: {water_count}")
    """
    # Define recognized water residue names (case-insensitive)
    water_residues = {'HOH', 'WAT'}

    # Initialize the count of water atoms
    water_atom_count = 0

    # Iterate through residues in the PDB file and count water atoms
    for residue in pdb.topology.residues():
        if residue.name.upper() in water_residues:
            water_atom_count += +1 

    return water_atom_count


pdb =PDBFile("/mnt/home/bturzo/ceph/Projects/Holography/03_Protein_Design/outputs/P_6_1_L_37_39_5J0K_4UOS/ForcedMonomer/MD/Mutation_9/R_1/energy_min.pdb")
print(get_water_count(pdb))


def read_params(filename):
    """
    Reads a list of CHARMM parameter files from a specified input file and returns a `CharmmParameterSet` object.

    This function parses an input file containing paths to parameter files, ignoring comments and blank lines.
    The parameter files are then loaded and combined into a single `CharmmParameterSet` object.

    Parameters
    ----------
    filename : str
        Path to the input file that contains a list of CHARMM parameter file paths.
        Each line should contain a single file path, with comments optionally following
        the `!` character.

    Returns
    -------
    params : CharmmParameterSet
        A `CharmmParameterSet` object containing the parameters from the listed files.

    Raises
    ------
    FileNotFoundError
        If the specified file does not exist.
    IOError
        If the file cannot be read due to permission or I/O errors.
    ValueError
        If no valid parameter files are specified in the input file.

    Notes
    -----
    - Lines starting with `!` or containing `!` after a file path are treated as comments and ignored.
    - Each valid, non-empty line is treated as a parameter file path.
    - `CharmmParameterSet` must be imported and available in the environment.

    Example
    -------
    >>> params = read_params('params_list.txt')
    >>> type(params)
    <class 'CharmmParameterSet'>

    File content (`params_list.txt`):
    ---------------------------------
    par_all36_prot.prm  ! Protein parameters
    par_all36_lipid.prm ! Lipid parameters
    par_all36_na.prm    ! Nucleic acid parameters

    The resulting `params` object will contain parameters from all listed files.

    """
    import os

    # Check if the input file exists
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"The file '{filename}' was not found.")

    parFiles = ()

    try:
        # Open and read the parameter file list
        with open(filename, 'r') as f:
            for line in f:
                print(f"Debug Line 1: {line}")
                # Remove comments after '!' and strip whitespace
                if '!' in line:
                    line = line.split('!')[0]
                parfile = line.strip()
                #print(f"Debug Line 1: {parfile}")
                # Add non-empty lines to the parameter file list
                if len(parfile) != 0:
                    parFiles += (parfile,)
                    #print(f"Debug Line 2: {parfile}")
        # Check if at least one valid parameter file was specified
        if not parFiles:
            raise ValueError(f"No valid parameter files were specified in '{filename}'.")

        # Load the parameter files using CharmmParameterSet
        print(f"Debug Line 3: {parFiles}")
        params = CharmmParameterSet(*parFiles)

    except IOError as e:
        raise IOError(f"Error reading file '{filename}': {e}")

    return params


print(read_params("./toppar.str"))
