# -*- coding: utf-8 -*-
"""
@author: Bargeen
"""
import sys
import argparse
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
from math import *
import os

parser = argparse.ArgumentParser()
parser.add_argument('-ipdb', dest='input_pdbfile', type=str, help='Input coordinate pdb file')
parser.add_argument('-itoppar', dest='path_to_topparstr', type=str, help='Input coordinate pdb file')
parser.add_argument('-ipsf', dest='path_to_psffile', type=str, help='Input coordinate pdb file')
parser.add_argument('-icrdf', dest='path_to_crdfile', type=str, help='Input coordinate pdb file')


parser.add_argument('-restart', dest='restart', help='Flag to restart simulation', action='store_true', default=False)
parser.add_argument('-solvent', dest='add_solvent', help='Flag to use solvent', action='store_true', default=False)
parser.add_argument('-charmmfff', dest='charmm_files_provided', help='If charmm files are provide use this flag to mention that', action='store_true', default=False)
parser.add_argument('-irstno',dest='currnt_restart_number', type=int, default=0, help='Restart number, by default it is set to zero')
parser.add_argument('-eqsteps',dest='equil_steps', type=int, default=50000, help='Number of Equilibration Steps')
parser.add_argument('-estore_report',dest='equil_store_report_interval', type=int, default=1000, help='At what steps of interval to store report for equilibration')
parser.add_argument('-pstore_report',dest='store_report_interval', type=int, default=5000, help='At what steps of interval to store report for production')
parser.add_argument('-psteps',dest='production_steps', type=int, default=125000000, help='Number of Production steps')
parser.add_argument('-simtep',dest='simulation_temperature', type=int, default=300, help='MD Simulation temperature')
parser.add_argument('-sim_type', dest='sim_run', type=str, help="If job is restarting, is it nvt or npt? if none is provided after restart flag, sys will exit.")
args = parser.parse_args()
# Args def
input_pdbfile               = args.input_pdbfile
restart                     = args.restart
add_solvent                 = args.add_solvent
currnt_restart_number       = args.currnt_restart_number
equil_steps                 = args.equil_steps
equil_store_report_interval = args.equil_store_report_interval
store_report_interval       = args.store_report_interval
production_steps            = args.production_steps
simulation_temperature      = args.simulation_temperature
sim_run                     = args.sim_run
charmm_files_provided       = args.charmm_files_provided
path_to_topparstr           = args.path_to_topparstr
path_to_psffile             = args.path_to_psffile
path_to_crdfile             = args.path_to_crdfile


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

    # Initialize the count of water residues
    water_count = 0

    # Iterate through residues in the PDB file and count water atoms
    for residue in pdb.topology.residues():
        if residue.name.upper() in water_residues:
            water_count += +1

    return water_count


def read_top(filename, fftype='CHARMM'):
    """
    Read a topology file and return the appropriate topology object 
    based on the specified force field type.

    Parameters:
    ----------
    filename : str
        The path to the topology file to be read.
    
    fftype : str, optional
        The type of force field used to generate the topology file.
        Supported options are:
        - 'CHARMM' : Reads a CHARMM PSF file using CharmmPsfFile.
        - 'AMBER'  : Reads an AMBER topology file using AmberPrmtopFile.
        Default is 'CHARMM'.

    Returns:
    -------
    top : CharmmPsfFile or AmberPrmtopFile
        An object representing the topology file, either as a 
        CharmmPsfFile or AmberPrmtopFile depending on the force field type.

    Raises:
    ------
    ValueError
        If an unsupported force field type is provided.

    Examples:
    --------
    >>> top = read_top('protein.psf', fftype='CHARMM')
    >>> top = read_top('protein.prmtop', fftype='AMBER')
    """
    # Check if the file exists
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"File '{filename}' not found.")
    try:
        if   fftype == 'CHARMM': top = CharmmPsfFile(filename)
        elif fftype == 'AMBER':  top = AmberPrmtopFile(filename)
        else:
            raise ValueError(f"Unsupported force field type: {fftype}")
    except Exception as e:
        raise IOError(f"Error reading file '{filename}': {str(e)}")

    return top

def read_crd(filename, fftype='CHARMM'):
    """
    Reads a coordinate file and returns the corresponding object based on the specified force field type.

    Parameters:
    ----------
    filename : str
        Path to the coordinate file to be read.
    fftype : str, optional
        Force field type used to interpret the coordinate file. 
        Supported options:
        - 'CHARMM' (default): Reads a CHARMM coordinate file using CharmmCrdFile.
        - 'AMBER': Reads an AMBER coordinate file using AmberInpcrdFile.

    Returns:
    -------
    crd : CharmmCrdFile or AmberInpcrdFile
        Coordinate file object corresponding to the specified force field.

    Raises:
    ------
    FileNotFoundError
        If the specified file does not exist.
    ValueError
        If an unsupported force field type is provided.
    IOError
        If the file cannot be read due to permission or I/O errors.

    Example:
    -------
    crd = read_crd('input.crd', fftype='CHARMM')
    """
    # Check if file exists and is accessible
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"The file '{filename}' was not found.")
    try:
        if fftype == 'CHARMM': crd = CharmmCrdFile(filename)
        elif fftype == 'AMBER': crd = AmberInpcrdFile(filename)
        else:
            raise ValueError(f"Unsupported force field type: '{fftype}'. Choose either 'CHARMM' or 'AMBER'.")
    except IOError as e:
        raise IOError(f"Error reading file '{filename}': {e}")
    
    return crd

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
    # Check if the input file exists
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"The file '{filename}' was not found.")

    parFiles = ()

    try:
        # Open and read the parameter file list
        with open(filename, 'r') as f:
            for line in f:
                #print(f"Debug Line 1: {line}")
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
        #print(f"Debug Line 3: {parFiles}")
        params = CharmmParameterSet(*parFiles)

    except IOError as e:
        raise IOError(f"Error reading file '{filename}': {e}")

    return params

# This function obtained from Charmm-Guii is ignored
# Because here we are generating the PBC box ourselves
# In the next function gen_box
#def read_box(psf, filename):
#    try:
#        sysinfo = json.load(open(filename, 'r'))
#        boxlx, boxly, boxlz = map(float, sysinfo['dimensions'][:3])
#    except:
#        for line in open(filename, 'r'):
#            segments = line.split('=')
#            if segments[0].strip() == "BOXLX": boxlx = float(segments[1])
#            if segments[0].strip() == "BOXLY": boxly = float(segments[1])
#            if segments[0].strip() == "BOXLZ": boxlz = float(segments[1])
#    psf.setBox(boxlx*angstroms, boxly*angstroms, boxlz*angstroms)
#    return psf

def gen_box(psf, crd):
    coords = crd.positions

    min_crds = [coords[0][0], coords[0][1], coords[0][2]]
    max_crds = [coords[0][0], coords[0][1], coords[0][2]]

    for coord in coords:
        min_crds[0] = min(min_crds[0], coord[0])
        min_crds[1] = min(min_crds[1], coord[1])
        min_crds[2] = min(min_crds[2], coord[2])
        max_crds[0] = max(max_crds[0], coord[0])
        max_crds[1] = max(max_crds[1], coord[1])
        max_crds[2] = max(max_crds[2], coord[2])

    boxlx = max_crds[0]-min_crds[0]
    boxly = max_crds[1]-min_crds[1]
    boxlz = max_crds[2]-min_crds[2]

    psf.setBox(boxlx, boxly, boxlz)
    return psf

def gen_box(psf, crd):
    """
    Generates and sets box dimensions based on the minimum and maximum coordinates from a coordinate file.

    This function computes the bounding box dimensions of a molecular system by determining
    the minimum and maximum coordinates from the given coordinate data (`crd`), and then sets
    the box dimensions for the `psf` (Protein Structure File) object.

    Parameters
    ----------
    psf : CharmmPsfFile or similar
        The PSF object whose box dimensions will be set.
    crd : CharmmCrdFile, AmberInpcrdFile, or similar
        The coordinate file object containing atomic positions.
        The object should have a `positions` attribute that provides atomic coordinates.

    Returns
    -------
    psf : CharmmPsfFile or similar
        The PSF object with the computed box dimensions set.

    Raises
    ------
    AttributeError
        If the `crd` object does not have a `positions` attribute.
    ValueError
        If the coordinate data is empty or the dimensions cannot be computed.

    Notes
    -----
    - The box dimensions are calculated by determining the minimum and maximum
      x, y, and z coordinates across all atoms in the system.
    - `psf.setBox()` should be a valid method to apply the box dimensions.
    - The units of the box dimensions are assumed to be consistent with the input coordinates.

    Example
    -------
    psf = CharmmPsfFile('input.psf')
    crd = CharmmCrdFile('input.crd')
    psf = gen_box(psf, crd)
    print(psf.boxLengths)  
    """
    # Check if the crd object has a 'positions' attribute
    if not hasattr(crd, 'positions'):
        raise AttributeError("The coordinate object does not have a 'positions' attribute.")

    coords = crd.positions

    # Check if the coordinate list is non-empty
    if len(coords) == 0:
        raise ValueError("The coordinate data is empty. Cannot generate box dimensions.")

    # Initialize min and max coordinates with the first atom's coordinates
    min_crds = [coords[0][0], coords[0][1], coords[0][2]]
    max_crds = [coords[0][0], coords[0][1], coords[0][2]]

    # Iterate through all coordinates to compute min and max bounds
    for coord in coords:
        min_crds[0] = min(min_crds[0], coord[0])
        min_crds[1] = min(min_crds[1], coord[1])
        min_crds[2] = min(min_crds[2], coord[2])
        max_crds[0] = max(max_crds[0], coord[0])
        max_crds[1] = max(max_crds[1], coord[1])
        max_crds[2] = max(max_crds[2], coord[2])

    # Calculate the box dimensions
    boxlx = max_crds[0] - min_crds[0]
    boxly = max_crds[1] - min_crds[1]
    boxlz = max_crds[2] - min_crds[2]

    # Set the box dimensions in the PSF object
    psf.setBox(boxlx, boxly, boxlz)
    return psf

def rewrap(simulation):
    """
    Rewraps molecular coordinates in a periodic simulation box to ensure bonded atoms remain close together.

    This function analyzes the positions of atoms in a periodic box, detects bonds that cross the box boundaries,
    and translates the atoms in the residue of the second bonded atom (`res2`) to keep it close to the first atom (`res1`).
    This process ensures that molecules remain wrapped properly inside the periodic simulation box.

    Parameters
    ----------
    simulation :
        The OpenMM simulation object that contains:
        - The topology with bond information.
        - The context that holds the current state and positions of the atoms.
        - The periodic box vectors defining the simulation box dimensions.

    Returns
    -------
    simulation : 
        The updated simulation object with rewrapped atom positions.

    Raises
    ------
    ValueError
        If no bonds are found in the topology or positions cannot be retrieved properly.

    Notes
    -----
    - Periodic boundary conditions can lead to atoms in bonded residues appearing on opposite sides of the box.
    - This function rewraps such atoms to keep bonded residues close together.
    - It calculates the center of the system and uses box dimensions to correct atomic positions.

    Example
    -------
    simulation = Simulation(topology, system, integrator)
    simulation.context.setPositions(initial_positions)
    simulation = rewrap(simulation)
    """
    bonds = simulation.topology.bonds()
    positions = simulation.context.getState(getPositions=True).getPositions()
    box = simulation.context.getState().getPeriodicBoxVectors()
    boxlx = box[0][0]/angstrom
    boxly = box[1][1]/angstrom
    boxlz = box[2][2]/angstrom
    min_crds = [positions[0][0]/angstrom, positions[0][1]/angstrom, positions[0][2]/angstrom]
    max_crds = [positions[0][0]/angstrom, positions[0][1]/angstrom, positions[0][2]/angstrom]
    for position in positions:
        min_crds[0] = min(min_crds[0], position[0]/angstrom)
        min_crds[1] = min(min_crds[1], position[1]/angstrom)
        min_crds[2] = min(min_crds[2], position[2]/angstrom)
        max_crds[0] = max(max_crds[0], position[0]/angstrom)
        max_crds[1] = max(max_crds[1], position[1]/angstrom)
        max_crds[2] = max(max_crds[2], position[2]/angstrom)
    xcen = (max_crds[0] + min_crds[0]) / 2.0
    ycen = (max_crds[1] + min_crds[1]) / 2.0
    zcen = (max_crds[2] + min_crds[2]) / 2.0
    for bond in bonds:
        atom1 = bond[0]
        atom2 = bond[1]
        atom1id = atom1.index
        atom2id = atom2.index
        res1 = atom1.residue
        res2 = atom2.residue
        x1, y1, z1 = positions[atom1id]
        x2, y2, z2 = positions[atom2id]
        dx = fabs(x1/angstrom - x2/angstrom)
        dy = fabs(y1/angstrom - y2/angstrom)
        dz = fabs(z1/angstrom - z2/angstrom)
        if dx > boxlx/2 or dy > boxly/2 or dz > boxlz/2:
            for atom in res2.atoms():
                oldx = positions[atom.index][0]/angstrom
                oldy = positions[atom.index][1]/angstrom
                oldz = positions[atom.index][2]/angstrom
                if dx > boxlx/2.0:
                    if oldx < xcen: newx = oldx + boxlx
                    else: newx = oldx - boxlx
                else:
                    newx = oldx
                if dy > boxly/2.0:
                    if oldy < ycen: newy = oldy + boxly
                    else: newy = oldy - boxly
                else:
                    newy = oldy
                if dz > boxlz/2.0:
                    if oldz < zcen: newz = oldz + boxlz
                    else: newz = oldz - boxlz
                else:
                    newz = oldz
                new_position = Vec3(newx, newy, newz)
                positions[atom.index] = Quantity(new_position, angstroms)
    simulation.context.setPositions(positions)
    return simulation

def set_compute_system(platform):
    """
    Selects and configures the appropriate OpenMM compute platform based on user input or available platforms.

    This function checks the available platforms in the OpenMM installation and either uses the
    specified platform or defaults to the best available platform in a predefined order (`CUDA`, `OpenCL`, `CPU`).
    It also configures platform properties (such as `CudaPrecision`) if applicable.

    Parameters
    ----------
    platform : str or None
        The name of the desired OpenMM platform.
        - If `platform` is provided, the function tries to use that platform.
        - If `platform` is `None`, the function selects the first available platform from the
          list `['CUDA', 'OpenCL', 'CPU']` in that order.

    Returns
    -------
    list
        A list containing:
        - platform (Platform) : The selected OpenMM platform.
        - prop (dict) : A dictionary of platform-specific properties.
          - For CUDA, the property `CudaPrecision` is set to `single`.
          - For other platforms, an empty dictionary is returned.

    Raises
    ------
    SystemExit
        If the specified platform is not available or no valid platform is found.

    Notes
    -----
    - CUDA is preferred over OpenCL and CPU if available, providing better performance.
    - If `platform` is specified but not found, the program exits with an error message.
    - If no valid platform is found when `platform` is `None`, the program exits with an error message.

    Example
    -------
    platform, prop = set_compute_system('CUDA')
    Using platform: CUDA

    platform, prop = set_compute_system(None)
    Using platform: CUDA
    """

    DEFAULT_PLATFORMS = 'CUDA', 'OpenCL', 'CPU'
    # Get the list of available platforms
    enabled_platforms = [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]
    print(enabled_platforms)

    # If a platform is specified, try to use it
    if platform:
        if platform not in enabled_platforms:
            print("Unable to find OpenMM platform '{}'; exiting".format(platform), file=sys.stderr)
            sys.exit(1)
        platform = Platform.getPlatformByName(platform)
    else:
        # Attempt to use the best available platform in the preferred order
        for platform in DEFAULT_PLATFORMS:
            if platform in enabled_platforms:
                platform = Platform.getPlatformByName(platform)
                break
        # If no valid platform was found, exit
        if isinstance(platform, str):
            print("Unable to find any OpenMM platform; exiting", file=sys.stderr)
            sys.exit(1)

    print("Using platform:", platform.getName())
    
    # Set platform-specific properties if needed
    prop = dict(CudaPrecision='single') if platform.getName() == 'CUDA' else dict()
    return [platform, prop]

def save_simulation_results(simulation, output_pdb, state_xml, system_xml):
    """
    Saves the results of an OpenMM simulation, including atom positions, system state, and system configuration.

    This function extracts the current state of the simulation and writes the results to:
    - A PDB file containing atom positions.
    - An XML file containing the state information.
    - An XML file containing the system definition.

    Parameters
    ----------
    simulation : 
        The OpenMM simulation object containing the topology, system, and current state.
    output_pdb : str or Path
        Path to the output PDB file where the atomic positions will be saved.
    state_xml : str or Path
        Path to the output XML file where the simulation state (positions, velocities, and forces) will be saved.
    system_xml : str or Path
        Path to the output XML file where the system configuration (forces, particles, constraints, etc.) will be saved.

    Returns
    -------
    None
        The function writes the simulation results to the specified files and does not return any value.

    Notes
    -----
    - The simulation state includes information such as positions, velocities, forces, energy, and parameter values.
    - The resulting PDB file is human-readable and can be visualized with molecular visualization tools.
    - The XML files provide detailed descriptions of the system and state, which can be used to restart or analyze the simulation.

    Example
    -------
    simulation = Simulation(topology, system, integrator)
    simulation.step(1000)  # Run 1000 steps of the simulation
    save_simulation_results(simulation, 'output.pdb', 'state.xml', 'system.xml')
    Saving results...
    """

    print("saving results")

    # Extract the current state with positions, velocities, forces, and energies
    state = simulation.context.getState(
        getPositions=True,
        getVelocities=True,
        getForces=True,
        getEnergy=True,
        getParameters=True,
        enforcePeriodicBox=True
    )

    # Save positions to a PDB file
    with open(str(output_pdb), 'w') as f:
        PDBFile.writeFile(simulation.topology, state.getPositions(), f)

    # Save the state to an XML file
    with open(str(state_xml), 'w') as f:
        f.write(XmlSerializer.serialize(state))

    # Save the system to an XML file
    with open(str(system_xml), 'w') as f:
        f.write(XmlSerializer.serialize(simulation.system))


def read_simulation_result(state_xml,system_xml):
    """
    Reads and deserializes the system and state information from XML files.

    This function loads a previously saved OpenMM system and state from XML files,
    allowing for the recreation of a simulation's configuration and state.

    Parameters
    ----------
    state_xml : str or Path
        Path to the XML file containing the serialized state information.
        The state file typically contains details about positions, velocities, forces,
        and system parameters at a specific simulation step.
    system_xml : str or Path
        Path to the XML file containing the serialized system information.
        The system file includes details such as particles, forces, constraints, and 
        other system definitions.

    Returns
    -------
    list
        A list containing:
        - system : openmm.System
            The deserialized OpenMM system.
        - state : openmm.State
            The deserialized OpenMM state.

    Raises
    ------
    FileNotFoundError
        If either the `state_xml` or `system_xml` file is not found.
    IOError
        If the files cannot be read due to permission or other I/O issues.
    ValueError
        If the XML files are corrupted or cannot be deserialized.

    Notes
    -----
    - Deserialization of the XML files restores the system and state objects,
      which can be used to restart a simulation from a specific state.
    - `XmlSerializer.deserialize()` is used to read and deserialize the XML content.
    
    Example
    -------
    system, state = read_simulation_result('state.xml', 'system.xml')
    print(system.getNumParticles())
    output: 1000
    print(state.getTime())
    output: 50.0 ps
    """    
    with open(str(system_xml), 'r') as f:
        system = XmlSerializer.deserialize(f.read())
    with open(str(state_xml), 'r') as f:
        state = XmlSerializer.deserialize(f.read())
    return [system, state]

def energy_minimization(modeller, system, platform, prop):
    """
    Performs energy minimization on a molecular system using OpenMM.

    This function sets up a simulation with a Langevin integrator and performs energy
    minimization to relax the system and remove any steric clashes.

    Parameters
    ----------
    modeller : 
        The Modeller object containing the system topology and initial atomic positions.
    system : 
        The OpenMM System object defining the forces and constraints of the system.
    platform : 
        The OpenMM platform used to run the simulation (e.g., CUDA, OpenCL, CPU).
    prop : dict
        Platform-specific properties, such as `CudaPrecision` for CUDA.

    Returns
    -------
    simulation : openmm.app.Simulation
        The Simulation object after energy minimization, with updated atomic positions.

    Notes
    -----
    - The integrator used is `LangevinMiddleIntegrator` with:
      - Temperature: 300 K
      - Friction coefficient: 1/picosecond
      - Time step: 2 fs (0.002 picoseconds)
    - `simulation.minimizeEnergy()` uses OpenMM's default energy minimization method.
    - The resulting `simulation` can be further used to perform dynamics or save results.

    Example
    -------
    modeller = Modeller(topology, positions)
    system = forcefield.createSystem(modeller.topology)
    simulation = energy_minimization(modeller, system, platform, {'CudaPrecision': 'single'})
    state = simulation.context.getState(getPositions=True)
    minimized_positions = state.getPositions(asNumpy=True)
    """    
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    simulation = Simulation(modeller.topology, system, integrator, platform, prop)
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy()
    return simulation

def charmm_energy_minimization(psffile, crdfile, system, platform, prop):
    """
    Performs energy minimization for a CHARMM system using OpenMM.

    This function sets up a simulation for a system defined by a CHARMM PSF file and
    a coordinate (CRD) file, and performs energy minimization to eliminate steric clashes
    and bring the system to a local energy minimum.

    Parameters
    ----------
    psffile : openmm.app.CharmmPsfFile
        The CHARMM PSF (Protein Structure File) containing the system topology.
    crdfile : openmm.app.CharmmCrdFile
        The CHARMM CRD (Coordinate File) containing initial atomic positions.
    system : openmm.System
        The OpenMM System object that defines the forces and constraints of the system.
    platform : openmm.Platform
        The OpenMM platform used to run the simulation (e.g., CUDA, OpenCL, CPU).
    prop : dict
        Platform-specific properties, such as `CudaPrecision` for CUDA.

    Returns
    -------
    simulation : openmm.app.Simulation
        The Simulation object after energy minimization, with updated atomic positions.

    Notes
    -----
    - The integrator used is `LangevinMiddleIntegrator` with:
      - Temperature: 300 K
      - Friction coefficient: 1/picosecond
      - Time step: 2 fs (0.002 picoseconds)
    - `simulation.minimizeEnergy()` minimizes the system energy using OpenMM's default
      energy minimization algorithm.
    - The resulting `simulation` object can be used for further simulations or analysis.

    Example
    -------
    from openmm.app import CharmmPsfFile, CharmmCrdFile
    psf = CharmmPsfFile('input.psf')
    crd = CharmmCrdFile('input.crd')
    system = forcefield.createSystem(psf.topology)
    simulation = charmm_energy_minimization(psf, crd, system, platform, {'CudaPrecision': 'single'})
    state = simulation.context.getState(getPositions=True)
    minimized_positions = state.getPositions(asNumpy=True)
    """
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    simulation = Simulation(psffile.topology, system, integrator, platform, prop)
    simulation.context.setPositions(crdfile.positions)
    simulation.minimizeEnergy()
    return simulation


def sim_context(pdb, system, state, store_report_interval, nsteps, dcd_filename, platform, prop, append_bool=False):
    """
    Runs a molecular dynamics (MD) simulation in OpenMM using a given initial state.

    This function initializes a simulation context using a provided PDB structure, system definition,
    and state. It sets up a Langevin integrator and appends reporters to monitor and save simulation
    data, including a DCD trajectory, state information, and checkpoint files. The simulation is run
    for a specified number of steps.

    Parameters
    ----------
    pdb : PDBFile
        The PDB file containing the initial atomic positions and topology information.
    system : System
        The OpenMM System object defining the forces, constraints, and particles of the system.
    state : State
        The initial state of the system, containing positions, velocities, and box dimensions.
    store_report_interval : int
        Interval (in steps) at which simulation data is written to output files.
    nsteps : int
        The total number of simulation steps to be performed.
    dcd_filename : str
        Name of the output DCD file where atomic positions will be stored.
    platform : Platform
        The OpenMM platform used to run the simulation (e.g., CUDA, OpenCL, CPU).
    prop : dict
        Platform-specific properties, such as `CudaPrecision` for CUDA.
    append_bool : bool, optional
        Whether to append data to an existing DCD file. Defaults to `False`.

    Returns
    -------
    simulation : openmm.app.Simulation
        The Simulation object after running the specified number of steps, with updated atomic positions.

    Notes
    -----
    - The integrator used is `LangevinMiddleIntegrator` with:
      - Temperature: 300 K
      - Friction coefficient: 1/picosecond
      - Time step: 2 fs (0.002 picoseconds)
    - The simulation writes:
        - Trajectories to a DCD file (`{dcd_filename}.dcd`).
        - Simulation state information (energy and temperature) to `stdout` at regular intervals.
        - Checkpoints to `checkpoint.chk` every 200 steps.
    - The resulting `simulation` object can be used to extract state information or extend the simulation.

    Example
    -------
    from openmm.app import PDBFile
    pdb = PDBFile('input.pdb')
    system = forcefield.createSystem(pdb.topology)
    simulation = sim_context(pdb, system, state, 1000, 50000, 'output', platform, {'CudaPrecision': 'single'})
    state = simulation.context.getState(getPositions=True)
    positions = state.getPositions(asNumpy=True)
    """

    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    simulation = Simulation(pdb.topology, system, integrator, platform, prop)
    simulation.context.setState(state)
    simulation.reporters.append(DCDReporter(f'{dcd_filename}.dcd', int(store_report_interval), append=append_bool))
    simulation.reporters.append(StateDataReporter(stdout, int(store_report_interval), step=True, potentialEnergy=True, temperature=True))
    checkpointReporter = CheckpointReporter('checkpoint.chk', 200)
    simulation.reporters.append(checkpointReporter)
    simulation.step(int(nsteps))
    return simulation

def nvt_run(pdb, system, state, store_report_interval, nsteps, dcd_filename, platform, prop, append_bool):
    """
    Runs an NVT (constant number of particles, volume, and temperature) simulation using OpenMM.

    This function serves as a wrapper for `sim_context` to run a molecular dynamics simulation
    under NVT conditions with a given system, initial state, and specified simulation parameters.

    Parameters
    ----------
    pdb : openmm.app.PDBFile
        The PDB file containing the initial atomic positions and topology information.
    system : openmm.System
        The OpenMM System object defining the forces, constraints, and particles of the system.
    state : openmm.State
        The initial state of the system, containing positions, velocities, and box dimensions.
    store_report_interval : int
        Interval (in steps) at which simulation data is written to output files.
    nsteps : int
        The total number of simulation steps to be performed.
    dcd_filename : str
        Name of the output DCD file where atomic positions will be stored.
    platform : openmm.Platform
        The OpenMM platform used to run the simulation (e.g., CUDA, OpenCL, CPU).
    prop : dict
        Platform-specific properties, such as `CudaPrecision` for CUDA.
    append_bool : bool
        Whether to append data to an existing DCD file.

    Returns
    -------
    simulation : openmm.app.Simulation
        The Simulation object after running the specified number of steps, with updated atomic positions.

    Notes
    -----
    - This function uses `sim_context` to set up and run the simulation.
    - NVT conditions are maintained by using a Langevin integrator with:
      - Constant temperature (300 K).
      - Constant volume (no barostat applied).
    - The resulting `simulation` object can be used for further analysis or continuation.

    Example
    -------
    from openmm.app import PDBFile
    pdb = PDBFile('input.pdb')
    system = forcefield.createSystem(pdb.topology)
    simulation = nvt_run(pdb, system, state, 1000, 50000, 'output', platform, {'CudaPrecision': 'single'}, False)
    state = simulation.context.getState(getPositions=True)
    positions = state.getPositions(asNumpy=True)
    """
    return sim_context(pdb, system, state, store_report_interval, nsteps, dcd_filename, platform, prop, append_bool)

def npt_run(pdb, system, state, store_report_interval, nsteps, dcd_filename, platform, prop, append_bool):
    """
    Runs an NPT (constant number of particles, pressure, and temperature) simulation using OpenMM.

    This function modifies the given system to include a Monte Carlo barostat to maintain constant
    pressure and temperature, and then calls `sim_context` to run the simulation under NPT conditions.

    Parameters
    ----------
    pdb : openmm.app.PDBFile
        The PDB file containing the initial atomic positions and topology information.
    system : openmm.System
        The OpenMM System object defining the forces, constraints, and particles of the system.
    state : openmm.State
        The initial state of the system, containing positions, velocities, and box dimensions.
    store_report_interval : int
        Interval (in steps) at which simulation data is written to output files.
    nsteps : int
        The total number of simulation steps to be performed.
    dcd_filename : str
        Name of the output DCD file where atomic positions will be stored.
    platform : openmm.Platform
        The OpenMM platform used to run the simulation (e.g., CUDA, OpenCL, CPU).
    prop : dict
        Platform-specific properties, such as `CudaPrecision` for CUDA.
    append_bool : bool
        Whether to append data to an existing DCD file.

    Returns
    -------
    simulation : openmm.app.Simulation
        The Simulation object after running the specified number of steps, with updated atomic positions.

    Notes
    -----
    - A `MonteCarloBarostat` is added to the system to maintain constant pressure and temperature:
      - Pressure: 1 bar
      - Temperature: 300 K
    - This function uses `sim_context` to set up and run the simulation.
    - The resulting `simulation` object can be used for further analysis or continuation.

    Example
    -------
    from openmm.app import PDBFile
    pdb = PDBFile('input.pdb')
    system = forcefield.createSystem(pdb.topology)
    simulation = npt_run(pdb, system, state, 1000, 50000, 'output', platform, {'CudaPrecision': 'single'}, False)
    state = simulation.context.getState(getPositions=True)
    positions = state.getPositions(asNumpy=True)
    """
    system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))
    return sim_context(pdb, system, state, store_report_interval, nsteps, dcd_filename, platform, prop, append_bool)

def initialize_forcefield(solvent=True):
    """
    Initializes an AMBER14 force field for a molecular system with or without explicit solvent.

    This function loads the appropriate force field XML files from the AMBER14 package based on
    whether the system is solvated or not.

    Parameters
    ----------
    solvent : bool, optional
        If `True` (default), the system includes explicit solvent and the TIP3P water model is loaded.
        If `False`, only the protein force field is used without solvent.

    Returns
    -------
    forcefield : openmm.app.ForceField
        The initialized ForceField object containing the relevant force field parameters.

    Notes
    -----
    - The `amber14/protein.ff14SB.xml` file defines the force field parameters for proteins.
    - The `amber14/tip3p.xml` file defines the TIP3P water model for explicit solvation.
    - When `solvent` is set to `False`, only the protein force field is applied.

    Example
    -------
    ff = initialize_forcefield(solvent=True)
    print(ff.getGenerators())
    """    
    if solvent:
        return ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3p.xml')
    return ForceField('amber14/protein.ff14SB.xml')

def create_system(modeller, forcefield, solvent=True):
    """
    Creates an OpenMM System from a Modeller object and a specified force field.

    This function generates a System object that defines the interactions and constraints
    for a molecular system. It applies different settings for nonbonded interactions based on
    whether the system is solvated or not.

    Parameters
    ----------
    modeller : openmm.app.Modeller
        The Modeller object containing the system topology and atomic positions.
    forcefield : openmm.app.ForceField
        The ForceField object containing the parameters for the system.
    solvent : bool, optional
        If `True` (default), the system is treated as solvated and uses PME for nonbonded
        interactions with periodic boundary conditions.
        If `False`, the system is treated as non-solvated, using a cutoff for nonbonded
        interactions without periodic boundary conditions.

    Returns
    -------
    system : openmm.System
        The generated System object that defines the particles, forces, and constraints
        for the molecular system.

    Notes
    -----
    - For solvated systems (`solvent=True`):
        - Nonbonded method: PME (Particle Mesh Ewald)
        - Nonbonded cutoff: 1 nm
        - Constraints: HBonds
    - For non-solvated systems (`solvent=False`):
        - Nonbonded method: CutoffNonPeriodic
        - Nonbonded cutoff: 1 nm
        - Constraints: HBonds

    Example
    -------
    modeller = Modeller(topology, positions)
    forcefield = ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3p.xml')
    system = create_system(modeller, forcefield, solvent=True)
    print(system.getNumParticles())
    """
    if solvent:
        return forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
    return forcefield.createSystem(modeller.topology, nonbondedMethod=CutoffNonPeriodic, nonbondedCutoff=1*nanometer, constraints=HBonds)

def run_simulation(pdb, system, state, run_type, remaining_steps, output_prefix, platform, prop, restart):
    """
    Runs a molecular dynamics (MD) simulation in either NVT or NPT ensemble using OpenMM.

    This function selects and runs the appropriate simulation type (`nvt` or `npt`)
    based on the `run_type` parameter. It initializes and runs the simulation for a specified
    number of steps, saving the results to output files with a specified prefix.

    Parameters
    ----------
    pdb : openmm.app.PDBFile
        The PDB file containing the initial atomic positions and topology information.
    system : openmm.System
        The OpenMM System object defining the forces, constraints, and particles of the system.
    state : openmm.State
        The initial state of the system, containing positions, velocities, and box dimensions.
    run_type : str
        Type of simulation to run. Accepted values:
        - "nvt" : Constant Number of Particles, Volume, and Temperature.
        - "npt" : Constant Number of Particles, Pressure, and Temperature.
    remaining_steps : int
        The total number of simulation steps to be performed.
    output_prefix : str
        Prefix used to name the output files (DCD trajectory and state information).
    platform : openmm.Platform
        The OpenMM platform used to run the simulation (e.g., CUDA, OpenCL, CPU).
    prop : dict
        Platform-specific properties, such as `CudaPrecision` for CUDA.
    restart : bool
        Whether to append simulation data to an existing DCD file if continuing a previous simulation.

    Returns
    -------
    simulation : Simulation
        The Simulation object after running the specified number of steps, with updated atomic positions.

    Raises
    ------
    ValueError
        If `run_type` is not one of the accepted values (`"nvt"` or `"npt"`).

    Notes
    -----
    - `nvt_run()` runs the simulation in an NVT ensemble.
    - `npt_run()` runs the simulation in an NPT ensemble.
    - The resulting `simulation` object can be used to extract state information or continue the simulation.

    Example
    -------
    pdb = PDBFile('input.pdb')
    system = forcefield.createSystem(pdb.topology)
    simulation = run_simulation(pdb, system, state, 'nvt', 50000, 'output', platform, {'CudaPrecision': 'single'}, False)
    state = simulation.context.getState(getPositions=True)
    """
    if run_type == "nvt":
        return nvt_run(pdb, system, state, equil_store_report_interval, remaining_steps, output_prefix, platform, prop, restart)
    elif run_type == "npt":
        return npt_run(pdb, system, state, store_report_interval, remaining_steps, output_prefix, platform, prop, restart)
    else:
        raise ValueError("Invalid simulation run type. Must be 'nvt' or 'npt'.")

def restart_protocol(system, simulation, sim_run, equil_steps, production_steps, restart_number, platform, prop):
    state = simulation.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)
    with open("current_state.pdb", 'w') as pdb_f:
        PDBFile.writeFile(simulation.topology, state.getPositions(), pdb_f)
    pdb = PDBFile("current_state.pdb")
    current_step = state.getStepCount()
    remaining_steps = (equil_steps if sim_run == "nvt" else production_steps) - current_step
    print(f"Remaining steps for {sim_run}: {remaining_steps}")
    output_prefix = f"{sim_run}_equil_{restart_number}" if sim_run == "nvt" else f"npt_md_{restart_number}"
    simulation = run_simulation(pdb, system, state, sim_run, remaining_steps, output_prefix, platform, prop, True)
    save_simulation_results(simulation, f"{output_prefix}.pdb", f"state_{output_prefix}.xml", f"system_{output_prefix}.xml")


def charmm_protocol(psffile, crdfile, system, platform, prop, equil_steps, production_steps, restart_number):
    # Energy minimization
    simulation = charmm_energy_minimization(psffile, crdfile, system, platform, prop)
    simulation = rewrap(simulation)
    save_simulation_results(simulation, "energy_min.pdb", "state_min.xml", "system_min.xml")
    system, state = read_simulation_result("state_min.xml","system_min.xml")
    pdb = PDBFile('energy_min.pdb')
    # NVT phase
    output_prefix = f"nvt_equil_{restart_number}"
    simulation = nvt_run(pdb, system, state, equil_store_report_interval, equil_steps, output_prefix, platform, prop, False)
    simulation = rewrap(simulation)
    save_simulation_results(simulation, f"{output_prefix}.pdb", f"state_{output_prefix}.xml", f"system_{output_prefix}.xml")
    system, state = read_simulation_result(f"state_{output_prefix}.xml", f"system_{output_prefix}.xml")    
    pdb = PDBFile(f'{output_prefix}.pdb')
    # NPT phase
    output_prefix = f"npt_md_{restart_number}"
    simulation = npt_run(pdb, system, state, store_report_interval, production_steps, output_prefix, platform, prop, False)
    simulation = rewrap(simulation)
    save_simulation_results(simulation, f"{output_prefix}.pdb", f"state_{output_prefix}.xml", f"system_{output_prefix}.xml")

def amber_protocol(pdb_file, add_solvent, platform, prop, equil_steps, production_steps, restart_number):
    pdb = PDBFile(pdb_file)
    forcefield = initialize_forcefield(add_solvent)
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(forcefield)

    if add_solvent==False:
        print("Warning Running simulation in ~Vaccuum")

    if add_solvent:
        modeller.addSolvent(forcefield, padding=1.0*nanometers, ionicStrength=0.15*molar)

    system = create_system(modeller, forcefield, add_solvent)
    simulation = energy_minimization(modeller, system, platform, prop)

    save_simulation_results(simulation, "energy_min.pdb", "state_min.xml", "system_min.xml")
    system, state = read_simulation_result("state_min.xml", "system_min.xml")
    pdb = PDBFile('energy_min.pdb')

    # NVT phase
    output_prefix = f"nvt_equil_{restart_number}"
    simulation = nvt_run(pdb, system, state, equil_store_report_interval, equil_steps, output_prefix, platform, prop, False)
    save_simulation_results(simulation, f"{output_prefix}.pdb", f"state_{output_prefix}.xml", f"system_{output_prefix}.xml")
    
    if add_solvent:
        system, state = read_simulation_result(f"state_{output_prefix}.xml", f"system_{output_prefix}.xml")
        # NPT phase
        output_prefix = f"npt_md_{restart_number}"
        simulation = npt_run(pdb, system, state, store_report_interval, production_steps, output_prefix, platform, prop, False)
        save_simulation_results(simulation, f"{output_prefix}.pdb", f"state_{output_prefix}.xml", f"system_{output_prefix}.xml")



# System Config
platform, prop = set_compute_system("CUDA")

# Main control flow
if restart:
    integrator = LangevinMiddleIntegrator(simulation_temperature*kelvin, 1/picosecond, 0.002*picoseconds)
    if charmm_files_provided:
        psf = read_top(path_to_psffile)
        crd = read_crd(path_to_crdfile)
        psf = gen_box(psf, crd)
        params = read_params(path_to_topparstr)
        system = psf.createSystem(params, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
        simulation = Simulation(psf.topology, system, integrator, platform, prop)
        #simulation.context.setPositions(crd.positions) #not needed
    else:
        forcefield = initialize_forcefield()
        pdb = PDBFile('nvt_equil_0.pdb') #assuming nvt is completed.
        #integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
        modeller = Modeller(pdb.topology, pdb.positions)
        system = create_system(modeller, forcefield)
        simulation = Simulation(pdb.topology, system, integrator, platform, prop)

    with open('checkpoint.chk', 'rb') as f:
        simulation.context.loadCheckpoint(f.read())
    
    restart_protocol(system, simulation, sim_run, equil_steps, production_steps, currnt_restart_number, platform, prop)
else:
    if charmm_files_provided:
        psf = read_top(path_to_psffile)
        crd = read_crd(path_to_crdfile)
        psf = gen_box(psf, crd)
        params = read_params(path_to_topparstr)
        system = psf.createSystem(params, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
        charmm_protocol(psf, crd, system, platform, prop, equil_steps, production_steps, currnt_restart_number)

    else:
        amber_protocol(input_pdbfile, add_solvent, platform, prop, equil_steps, production_steps, currnt_restart_number)

