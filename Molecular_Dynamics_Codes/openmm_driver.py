# -*- coding: utf-8 -*-
"""
@author: Bargeen
"""
import argparse
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

parser = argparse.ArgumentParser()
parser.add_argument('-ipdb', dest='input_pdbfile', help='Input coordinate pdb file')
parser.add_argument('-restart', dest='restart', help='Flag to restart simulation', action='store_true', default=False)
parser.add_argument('-solvent', dest='add_solvent', help='Flag to use solvent', action='store_true', default=False)
parser.add_argument('-irstno',dest='currnt_restart_number', type=int, default=0, help='Restart number, by default it is set to zero')
parser.add_argument('-eqsteps',dest='equil_steps', type=int, default=50000, help='Number of Equilibration Steps')
parser.add_argument('-estore_report',dest='equil_store_report_interval', type=int, default=1000, help='At what steps of interval to store report for equilibration')
parser.add_argument('-pstore_report',dest='store_report_interval', type=int, default=5000, help='At what steps of interval to store report for production')
parser.add_argument('-psteps',dest='production_steps', type=int, default=125000000, help='Number of Production steps')
parser.add_argument('-simtep',dest='simulation_temperature', type=int, default=300, help='MD Simulation temperature')
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
    integrator = VerletIntegrator(0.001*picoseconds)
    simulation = Simulation(modeller.topology, system, integrator, platform, prop)
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy()
    return simulation

def sim_context(pdb, system, state, store_report_interval, nsteps, dcd_filename, platform, prop):
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    simulation = Simulation(pdb.topology, system, integrator, platform, prop)
    simulation.context.setState(state)
    simulation.reporters.append(DCDReporter(f'{dcd_filename}.dcd', int(store_report_interval)))
    simulation.reporters.append(StateDataReporter(stdout, int(store_report_interval), step=True, potentialEnergy=True, temperature=True))
    simulation.step(int(nsteps))
    return simulation

def nvt_run(pdb, system, state, store_report_interval, nsteps, dcd_filename, platform, prop):
    return sim_context(pdb, system, state, store_report_interval, nsteps, dcd_filename, platform, prop)

def npt_run(pdb, system, state, store_report_interval, nsteps, dcd_filename, platform, prop):
    system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))
    return sim_context(pdb, system, state, store_report_interval, nsteps, dcd_filename, platform, prop)

# System Config
platform, prop = set_compute_system("CUDA")

# MD Protocol Control Flow
if restart:
    restart_number = currnt_restart_number +1
    if add_solvent:
        print(f"Restarting solution phase simulation of {input_pdbfile}")
        pdb = PDBFile(f'npt_md_{currnt_restart_number}.pdb')
        system, state = read_simulation_result(f"state_npt_md_{currnt_restart_number}.xml", f"system_npt_md_{currnt_restart_number}.xml")
        simulation = npt_run(pdb, system, state, store_report_interval, production_steps, f"npt_md_{restart_number}", platform, prop)
        save_npt_results = save_simulation_results(simulation, f"npt_md_{restart_number}.pdb", f"state_npt_md_{restart_number}.xml", f"system_npt_md_{restart_number}.xml")
    else:
        print(f"Restarting vaccum phase simulation of {input_pdbfile}")
        pdb = PDBFile(f"nvt_md_{currnt_restart_number}.pdb")
        system, state = read_simulation_result(f"state_nvt_md_{currnt_restart_number}.xml", f"system_nvt_md_{currnt_restart_number}.xml")
        simulation = nvt_run(pdb, system, state, store_report_interval, production_steps, f"nvt_md_{restart_number}", platform, prop)
        save_nvt_results = save_simulation_results(simulation, f"nvt_md_{restart_number}.pdb", f"state_nvt_md_{restart_number}.xml", f"system_nvt_md_{restart_number}.xml")
else:
    pdb = PDBFile(f'{input_pdbfile}')
    if add_solvent:
        print(f"Starting a new solution phase simulation of {input_pdbfile}")
        #EM
        forcefield = ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3p.xml')
        modeller = Modeller(pdb.topology, pdb.positions)
        #modeller.addHydrogens(forcefield)
        modeller.addSolvent(forcefield, padding=1.0*nanometers, ionicStrength=0.15*molar)
        system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
        simulation = energy_minimization(modeller, system, platform, prop)
        save_energy_min_results = save_simulation_results(simulation, "energy_min.pdb", "state_min.xml", "system_min.xml")
        system, state = read_simulation_result("state_min.xml","system_min.xml")
        pdb = PDBFile('energy_min.pdb')                                                                                                                                      
        #NVT
        simulation = nvt_run(pdb, system, state, equil_store_report_interval, equil_steps, f"nvt_equil_{currnt_restart_number}", platform, prop)
        save_nvt_equil_results = save_simulation_results(simulation, f"nvt_equil_{currnt_restart_number}.pdb", f"state_nvt_equil_{currnt_restart_number}.xml", f"system_nvt_equil_{currnt_restart_number}.xml")
        pdb = PDBFile(f'nvt_equil_{currnt_restart_number}.pdb')
        system, state = read_simulation_result( f"state_nvt_equil_{currnt_restart_number}.xml",f"system_nvt_equil_{currnt_restart_number}.xml")
        #NPT
        simulation = npt_run(pdb, system, state, store_report_interval, production_steps, f"npt_md_{currnt_restart_number}", platform, prop)
        save_npt_results = save_simulation_results(simulation, f"npt_md_{currnt_restart_number}.pdb", f"state_npt_md_{currnt_restart_number}.xml", f"system_npt_md_{currnt_restart_number}.xml")
    else:
        print(f"Restarting a new vaccum phase simulation of {input_pdbfile}")
        #EM
        forcefield = ForceField('amber14/protein.ff14SB.xml')
        modeller = Modeller(pdb.topology, pdb.positions)
        #modeller.addHydrogens(forcefield)
        system = forcefield.createSystem(modeller.topology, nonbondedMethod=CutoffNonPeriodic, nonbondedCutoff=1*nanometer, constraints=HBonds)
        simulation = energy_minimization(modeller, system, platform, prop)
        save_energy_min_results = save_simulation_results(simulation, "energy_min.pdb", "state_min.xml", "system_min.xml")
        system, state = read_simulation_result("state_min.xml","system_min.xml")
        pdb = PDBFile('energy_min.pdb')
        #NVT
        simulation = nvt_run(pdb, system, state, store_report_interval, production_steps, f"nvt_md_{currnt_restart_number}", platform, prop) 
        save_nvt_results = save_simulation_results(simulation, f"nvt_md_{currnt_restart_number}.pdb", f"state_nvt_md_{currnt_restart_number}.xml", f"system_nvt_md_{currnt_restart_number}.xml")
