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
    # Initialize the count of water atoms
    water_atom_count = 0
    # Identify and count water atoms
    for residue in pdb.topology.residues():
        if residue.name == 'HOH' or residue.name == 'WAT':  # Check for water residue names
            for atom in residue.atoms():
                water_atom_count += 1

    return water_atom_count


def read_top(filename, fftype='CHARMM'):
    if   fftype == 'CHARMM': top = CharmmPsfFile(filename)
    elif fftype == 'AMBER':  top = AmberPrmtopFile(filename)
    return top

def read_crd(filename, fftype='CHARMM'):
    if   fftype == 'CHARMM': crd = CharmmCrdFile(filename)
    elif fftype == 'AMBER':  crd = AmberInpcrdFile(filename)
    return crd

def read_params(filename):                                                                                                             
    parFiles = ()
    for line in open(filename, 'r'):
        if '!' in line: line = line.split('!')[0]
        parfile = line.strip()
        if len(parfile) != 0: parFiles += ( parfile, )

    params = CharmmParameterSet( *parFiles )
    return params

def read_box(psf, filename):
    try:
        sysinfo = json.load(open(filename, 'r'))
        boxlx, boxly, boxlz = map(float, sysinfo['dimensions'][:3])
    except:
        for line in open(filename, 'r'):
            segments = line.split('=')
            if segments[0].strip() == "BOXLX": boxlx = float(segments[1])
            if segments[0].strip() == "BOXLY": boxly = float(segments[1])
            if segments[0].strip() == "BOXLZ": boxlz = float(segments[1])
    psf.setBox(boxlx*angstroms, boxly*angstroms, boxlz*angstroms)
    return psf

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

def rewrap(simulation):
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
    DEFAULT_PLATFORMS = 'CUDA', 'OpenCL', 'CPU' 
    enabled_platforms = [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())] 
    print(enabled_platforms)
    if platform:
        if not platform in enabled_platforms:
            print("Unable to find OpenMM platform '{}'; exiting".format(platform), file=sys.stderr)
            sys.exit(1)
        platform = Platform.getPlatformByName(platform)
    else:
        for platform in DEFAULT_PLATFORMS:
            if platform in enabled_platforms:
                platform = Platform.getPlatformByName(platform)
                break
        if isinstance(platform, str):
            print("Unable to find any OpenMM platform; exiting".format(platform), file=sys.stderr)
            sys.exit(1)
    print("Using platform:", platform.getName())
    prop = dict(CudaPrecision='single') if platform.getName() == 'CUDA' else dict()
    return [platform,prop]

def save_simulation_results(simulation, output_pdb, state_xml, system_xml):
    print("saving results")
    state = simulation.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)
    with open(str(output_pdb), 'w') as f:
        PDBFile.writeFile(simulation.topology, state.getPositions(), f)
    with open(str(state_xml), 'w') as f:
        f.write(XmlSerializer.serialize(state))
    with open(str(system_xml), 'w') as f:
        f.write(XmlSerializer.serialize(simulation.system))

def read_simulation_result(state_xml,system_xml):
    with open(str(system_xml), 'r') as f:
        system = XmlSerializer.deserialize(f.read())
    with open(str(state_xml), 'r') as f:
        state = XmlSerializer.deserialize(f.read())
    return [system, state]

def energy_minimization(modeller, system, platform, prop):
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    simulation = Simulation(modeller.topology, system, integrator, platform, prop)
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy()
    return simulation

def charmm_energy_minimization(psffile, crdfile, system, platform, prop):
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    simulation = Simulation(psffile.topology, system, integrator, platform, prop)
    simulation.context.setPositions(crdfile.positions)
    simulation.minimizeEnergy()
    return simulation


def sim_context(pdb, system, state, store_report_interval, nsteps, dcd_filename, platform, prop, append_bool=False):
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
    return sim_context(pdb, system, state, store_report_interval, nsteps, dcd_filename, platform, prop, append_bool)

def npt_run(pdb, system, state, store_report_interval, nsteps, dcd_filename, platform, prop, append_bool):
    system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))
    return sim_context(pdb, system, state, store_report_interval, nsteps, dcd_filename, platform, prop, append_bool)

def initialize_forcefield(solvent=True):
    if solvent:
        return ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3p.xml')
    return ForceField('amber14/protein.ff14SB.xml')

def create_system(modeller, forcefield, solvent=True):
    if solvent:
        return forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
    return forcefield.createSystem(modeller.topology, nonbondedMethod=CutoffNonPeriodic, nonbondedCutoff=1*nanometer, constraints=HBonds)

def run_simulation(pdb, system, state, run_type, remaining_steps, output_prefix, platform, prop, restart):
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

