# -*- coding: utf-8 -*-
"""
@author: Bargeen
"""
import argparse
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
from math import *

from openmmtools import states, mcmc, multistate
from openmmtools.states import ThermodynamicState, SamplerState

parser = argparse.ArgumentParser()
parser.add_argument('-topparstr', dest='path_to_topparstr', help='Input Topppar Stream file (the file with .str extentsion)')
parser.add_argument('-path_to_psf', dest='path_to_psf', help='Input psf file')
parser.add_argument('-path_to_crdfile', dest='path_to_crdfile', help='Input coordinate file')
parser.add_argument('-run_nvt_sim', dest='run_nvt_sim', help='Flag to run NVT simulation', action='store_true', default=False)
parser.add_argument('-run_trex_sim', dest='run_trex_sim', help='Flag to run Temperature Replica exchange simulation', action='store_true', default=False)
parser.add_argument('-restart', dest='restart', help='Flag to restart simulation', action='store_true', default=False)
parser.add_argument('-irstno',dest='currnt_restart_number', type=int, default=0, help='Restart number, by default it is set to zero')
parser.add_argument('-eqsteps',dest='equil_steps', type=int, default=5000000, help='Number of Equilibration Steps')
parser.add_argument('-n_iter',dest='n_iter', type=int, default=50000, help='Total number of Replica Exchange iteration')
parser.add_argument('-simtep',dest='simulation_temperature', type=int, default=310, help='MD Simulation temperature')
parser.add_argument('-pstore_report',dest='store_report_interval', type=int, default=5000, help='At what steps of interval to store report for production')
args = parser.parse_args()

# Args def
path_to_topparstr           = args.path_to_topparstr
path_to_psf                 = args.path_to_psf
path_to_crdfile             = args.path_to_crdfile
n_iter                      = args.n_iter
run_nvt_sim                 = args.run_nvt_sim
run_trex_sim                = args.run_trex_sim
restart                     = args.restart
currnt_restart_number       = args.currnt_restart_number
equil_steps                 = args.equil_steps
store_report_interval       = args.store_report_interval
simulation_temperature      = args.simulation_temperature


# TREMD DICTIONARY BASED ON WATER COUNT:
temp_dict = {
        20469: [310.00, 311.13, 312.27, 313.41, 314.55, 315.70, 316.85, 318.00, 319.16, 320.32, 321.49, 322.66, 323.83, 325.00, 326.18, 327.37, 328.55, 329.74, 330.93, 332.13, 333.33, 334.54, 335.74, 336.96, 338.17, 339.39, 340.61, 341.84, 343.07, 344.30, 345.54, 346.78, 348.03, 349.27, 350.53, 351.78, 353.04, 354.31, 355.58, 356.85, 358.13, 359.41, 360.69, 361.98, 363.27, 364.56, 365.86, 367.17, 368.47, 369.79, 371.10, 372.42, 373.75, 375.07, 376.40, 377.74, 379.08, 380.43, 381.78, 383.13, 384.49, 385.00], 20550: [310.00, 311.13, 312.26, 313.40, 314.54, 315.69, 316.84, 317.99, 319.14, 320.30, 321.47, 322.63, 323.80, 324.97, 326.15, 327.33, 328.52, 329.70, 330.89, 332.09, 333.29, 334.49, 335.69, 336.90, 338.12, 339.33, 340.55, 341.78, 343.00, 344.24, 345.47, 346.71, 347.95, 349.20, 350.45, 351.70, 352.97, 354.23, 355.50, 356.78, 358.04, 359.32, 360.59, 361.87, 363.15, 364.44, 365.73, 367.03, 368.34, 369.65, 370.96, 372.28, 373.60, 374.93, 376.26, 377.59, 378.93, 380.27, 381.61, 382.96, 384.31, 385.00],
        20475: [310.00, 311.14, 312.27, 313.41, 314.56, 315.70, 316.85, 318.01, 319.16, 320.32, 321.49, 322.66, 323.83, 325.00, 326.18, 327.36, 328.55, 329.74, 330.93, 332.13, 333.33, 334.53, 335.74, 336.95, 338.17, 339.40, 340.62, 341.85, 343.08, 344.31, 345.55, 346.79, 348.04, 349.29, 350.54, 351.79, 353.05, 354.32, 355.59, 356.86, 358.14, 359.42, 360.70, 361.99, 363.28, 364.58, 365.87, 367.18, 368.49, 369.80, 371.11, 372.43, 373.76, 375.08, 376.42, 377.75, 379.09, 380.43, 381.78, 383.13, 384.49, 385.00], 19389: [310.00, 311.16, 312.33, 313.50, 314.68, 315.86, 317.04, 318.22, 319.41, 320.61, 321.80, 323.00, 324.21, 325.42, 326.63, 327.84, 329.06, 330.29, 331.51, 332.75, 333.98, 335.22, 336.46, 337.71, 338.96, 340.21, 341.47, 342.74, 344.00, 345.27, 346.55, 347.82, 349.11, 350.39, 351.68, 352.98, 354.28, 355.58, 356.89, 358.20, 359.51, 360.83, 362.15, 363.48, 364.81, 366.15, 367.49, 368.83, 370.18, 371.53, 372.89, 374.25, 375.62, 376.99, 378.36, 379.74, 381.12, 382.50, 383.89, 385.00],
        19470: [310.00, 311.16, 312.33, 313.49, 314.67, 315.84, 317.02, 318.20, 319.39, 320.58, 321.78, 322.97, 324.18, 325.38, 326.59, 327.80, 329.02, 330.24, 331.47, 332.70, 333.93, 335.17, 336.40, 337.65, 338.89, 340.15, 341.40, 342.66, 343.92, 345.19, 346.46, 347.74, 349.02, 350.30, 351.59, 352.88, 354.17, 355.47, 356.78, 358.08, 359.40, 360.71, 362.03, 363.36, 364.68, 366.02, 367.35, 368.69, 370.04, 371.39, 372.74, 374.10, 375.46, 376.83, 378.20, 379.51, 380.89, 382.28, 383.66, 385.00],21921: [310.00, 311.10, 312.19, 313.30, 314.40, 315.51, 316.62, 317.74, 318.85, 319.97, 321.10, 322.23, 323.36, 324.49, 325.63, 326.77, 327.92, 329.06, 330.22, 331.37, 332.53, 333.69, 334.86, 336.03, 337.20, 338.37, 339.55, 340.73, 341.92, 343.11, 344.30, 345.50, 346.70, 347.90, 349.11, 350.32, 351.53, 352.75, 353.97, 355.20, 356.43, 357.66, 358.89, 360.13, 361.38, 362.62, 363.87, 365.13, 366.39, 367.65, 368.91, 370.18, 371.45, 372.73, 374.01, 375.29, 376.58, 377.87, 379.17, 380.47, 381.77, 383.08, 384.39, 385.00],
        19416: [310.00, 311.16, 312.33, 313.50, 314.67, 315.85, 317.03, 318.22, 319.41, 320.60, 321.79, 322.99, 324.20, 325.41, 326.62, 327.83, 329.07, 330.29, 331.52, 332.74, 333.97, 335.21, 336.45, 337.70, 338.95, 340.20, 341.46, 342.72, 343.97, 345.24, 346.51, 347.79, 349.09, 350.37, 351.66, 352.96, 354.25, 355.56, 356.86, 358.17, 359.48, 360.80, 362.12, 363.45, 364.78, 366.12, 367.46, 368.80, 370.15, 371.50, 372.86, 374.22, 375.58, 376.95, 378.32, 379.70, 381.08, 382.47, 383.86, 385.00], 19401: [310.00, 311.16, 312.33, 313.50, 314.68, 315.85, 317.04, 318.22, 319.41, 320.60, 321.80, 323.00, 324.21, 325.41, 326.62, 327.84, 329.06, 330.28, 331.51, 332.74, 333.97, 335.21, 336.46, 337.70, 338.95, 340.21, 341.47, 342.73, 343.99, 345.26, 346.54, 347.82, 349.10, 350.38, 351.67, 352.97, 354.27, 355.57, 356.87, 358.18, 359.50, 360.82, 362.14, 363.47, 364.80, 366.13, 367.47, 368.80, 370.15, 371.50, 372.86, 374.22, 375.58, 376.96, 378.33, 379.71, 381.09, 382.47, 383.87, 385.00],
        20532: [310.00, 311.13, 312.27, 313.41, 314.55, 315.69, 316.84, 317.99, 319.15, 320.31, 321.47, 322.64, 323.80, 324.98, 326.16, 327.34, 328.52, 329.71, 330.90, 332.09, 333.29, 334.49, 335.70, 336.91, 338.12, 339.34, 340.56, 341.78, 343.01, 344.24, 345.48, 346.72, 347.96, 349.21, 350.46, 351.71, 352.97, 354.24, 355.50, 356.77, 358.05, 359.32, 360.61, 361.89, 363.18, 364.48, 365.77, 367.08, 368.38, 369.69, 371.00, 372.32, 373.64, 374.97, 376.30, 377.63, 378.97, 380.31, 381.66, 383.01, 384.37, 385.00], 21906: [310.00, 311.10, 312.20, 313.30, 314.40, 315.51, 316.62, 317.74, 318.86, 319.98, 321.17, 322.29, 323.43, 324.56, 325.70, 326.84, 327.99, 329.14, 330.29, 331.44, 332.60, 333.76, 334.93, 336.10, 337.27, 338.45, 339.62, 340.80, 341.99, 343.18, 344.38, 345.57, 346.77, 347.98, 349.18, 350.40, 351.61, 352.82, 354.04, 355.26, 356.49, 357.73, 358.96, 360.20, 361.45, 362.69, 363.95, 365.20, 366.46, 367.72, 368.99, 370.26, 371.53, 372.81, 374.11, 375.40, 376.68, 377.97, 379.27, 380.57, 381.87, 383.18, 384.49, 385.00],
        20484: [310.00, 311.13, 312.26, 313.40, 314.55, 315.69, 316.84, 318.00, 319.16, 320.32, 321.48, 322.65, 323.82, 325.00, 326.18, 327.36, 328.54, 329.73, 330.93, 332.12, 333.32, 334.53, 335.73, 336.95, 338.16, 339.37, 340.59, 341.82, 343.05, 344.28, 345.52, 346.76, 348.00, 349.25, 350.50, 351.76, 353.02, 354.28, 355.55, 356.82, 358.10, 359.38, 360.66, 361.95, 363.24, 364.54, 365.84, 367.14, 368.45, 369.76, 371.07, 372.39, 373.72, 375.04, 376.37, 377.71, 379.05, 380.38, 381.73, 383.08, 384.44, 385.00], 20535: [310.00, 311.13, 312.26, 313.40, 314.54, 315.69, 316.84, 317.99, 319.15, 320.31, 321.47, 322.64, 323.81, 324.98, 326.16, 327.34, 328.52, 329.74, 330.93, 332.13, 333.32, 334.53, 335.73, 336.94, 338.15, 339.37, 340.59, 341.82, 343.05, 344.28, 345.51, 346.75, 348.00, 349.24, 350.49, 351.75, 353.01, 354.27, 355.54, 356.81, 358.08, 359.36, 360.64, 361.93, 363.22, 364.51, 365.81, 367.11, 368.41, 369.72, 371.04, 372.35, 373.68, 375.00, 376.33, 377.66, 379.00, 380.34, 381.69, 383.04, 384.39, 385.00],
        20505: [310.00, 311.13, 312.27, 313.41, 314.55, 315.69, 316.84, 318.00, 319.17, 320.33, 321.49, 322.66, 323.83, 325.00, 326.18, 327.36, 328.55, 329.74, 330.93, 332.12, 333.32, 334.53, 335.74, 336.95, 338.16, 339.38, 340.60, 341.83, 343.06, 344.29, 345.53, 346.77, 348.01, 349.26, 350.51, 351.77, 353.03, 354.29, 355.56, 356.83, 358.10, 359.38, 360.67, 361.95, 363.24, 364.54, 365.84, 367.14, 368.44, 369.75, 371.07, 372.39, 373.71, 375.04, 376.37, 377.70, 379.04, 380.38, 381.73, 383.08, 384.44, 385.00], 20496: [310.00, 311.13, 312.27, 313.41, 314.55, 315.70, 316.85, 318.00, 319.16, 320.32, 321.48, 322.65, 323.82, 324.99, 326.17, 327.35, 328.54, 329.73, 330.92, 332.12, 333.32, 334.52, 335.73, 336.94, 338.15, 339.37, 340.59, 341.82, 343.05, 344.28, 345.52, 346.76, 348.00, 349.25, 350.50, 351.76, 353.02, 354.28, 355.55, 356.82, 358.09, 359.37, 360.65, 361.94, 363.23, 364.52, 365.82, 367.12, 368.43, 369.74, 371.06, 372.37, 373.70, 375.02, 376.36, 377.69, 379.03, 380.37, 381.72, 383.07, 384.42, 385.00],
        20469: [310.00, 311.13, 312.27, 313.41, 314.55, 315.70, 316.85, 318.00, 319.16, 320.32, 321.49, 322.66, 323.83, 325.00, 326.18, 327.37, 328.55, 329.74, 330.93, 332.13, 333.33, 334.54, 335.74, 336.96, 338.17, 339.39, 340.61, 341.84, 343.07, 344.30, 345.54, 346.78, 348.03, 349.27, 350.53, 351.78, 353.04, 354.31, 355.58, 356.85, 358.13, 359.41, 360.69, 361.98, 363.27, 364.56, 365.86, 367.17, 368.47, 369.79, 371.10, 372.42, 373.75, 375.07, 376.40, 377.74, 379.08, 380.43, 381.78, 383.13, 384.49, 385.00], 20493: [310.00, 311.13, 312.27, 313.41, 314.55, 315.70, 316.85, 318.00, 319.16, 320.32, 321.48, 322.65, 323.82, 325.00, 326.17, 327.36, 328.54, 329.73, 330.92, 332.12, 333.32, 334.52, 335.73, 336.94, 338.15, 339.37, 340.59, 341.82, 343.05, 344.28, 345.52, 346.76, 348.00, 349.25, 350.51, 351.76, 353.02, 354.29, 355.55, 356.82, 358.10, 359.38, 360.66, 361.95, 363.24, 364.53, 365.83, 367.13, 368.44, 369.75, 371.07, 372.39, 373.71, 375.04, 376.37, 377.70, 379.04, 380.38, 381.73, 383.08, 384.43, 385.00],
        21870: [310.00, 311.10, 312.20, 313.30, 314.41, 315.51, 316.63, 317.74, 318.86, 319.99, 321.11, 322.24, 323.37, 324.51, 325.65, 326.79, 327.93, 329.08, 330.24, 331.41, 332.57, 333.73, 334.90, 336.07, 337.24, 338.42, 339.60, 340.78, 342.00, 343.19, 344.39, 345.59, 346.79, 347.99, 349.20, 350.41, 351.63, 352.85, 354.07, 355.30, 356.53, 357.77, 359.00, 360.25, 361.49, 362.74, 363.99, 365.25, 366.51, 367.77, 369.04, 370.31, 371.59, 372.87, 374.15, 375.44, 376.73, 378.02, 379.32, 380.62, 381.92, 383.23, 384.55, 385.00], 20496: [310.00, 311.13, 312.27, 313.41, 314.55, 315.70, 316.85, 318.00, 319.16, 320.32, 321.48, 322.65, 323.82, 324.99, 326.17, 327.35, 328.54, 329.73, 330.92, 332.12, 333.32, 334.52, 335.73, 336.94, 338.15, 339.37, 340.59, 341.82, 343.05, 344.28, 345.52, 346.76, 348.00, 349.25, 350.50, 351.76, 353.02, 354.28, 355.55, 356.82, 358.09, 359.37, 360.65, 361.94, 363.23, 364.52, 365.82, 367.12, 368.43, 369.74, 371.06, 372.37, 373.70, 375.02, 376.36, 377.69, 379.03, 380.37, 381.72, 383.07, 384.42, 385.00],
        19437: [310.00, 311.16, 312.33, 313.50, 314.67, 315.85, 317.03, 318.21, 319.40, 320.59, 321.79, 322.99, 324.19, 325.39, 326.60, 327.82, 329.04, 330.26, 331.48, 332.71, 333.95, 335.19, 336.43, 337.67, 338.92, 340.17, 341.43, 342.69, 343.96, 345.22, 346.50, 347.77, 349.05, 350.34, 351.63, 352.92, 354.21, 355.51, 356.82, 358.13, 359.44, 360.76, 362.08, 363.40, 364.73, 366.07, 367.41, 368.75, 370.09, 371.44, 372.80, 374.15, 375.52, 376.88, 378.26, 379.63, 381.01, 382.40, 383.79, 385.00]}


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

def energy_minimization(psffile, crdfile, system, platform, prop):
    integrator = VerletIntegrator(0.001*picoseconds)
    simulation = Simulation(psffile.topology, system, integrator, platform, prop)
    simulation.context.setPositions(crdfile.positions)
    simulation.minimizeEnergy()
    return simulation

def sim_context(pdb, system, state, store_report_interval, nsteps, dcd_filename, platform, prop):
    integrator = LangevinMiddleIntegrator(simulation_temperature*kelvin, 1/picosecond, 0.002*picoseconds)
    simulation = Simulation(pdb.topology, system, integrator, platform, prop)
    simulation.context.setState(state)
    simulation.reporters.append(DCDReporter(f'{dcd_filename}.dcd', int(store_report_interval)))
    simulation.reporters.append(StateDataReporter(stdout, int(store_report_interval), step=True, potentialEnergy=True, temperature=True))
    simulation.step(int(nsteps))
    return simulation

def nvt_run(pdb, system, state, store_report_interval, nsteps, dcd_filename, platform, prop):
    return sim_context(pdb, system, state, store_report_interval, nsteps, dcd_filename, platform, prop)

def npt_run(pdb, system, state, store_report_interval, nsteps, dcd_filename, platform, prop):
    system.addForce(MonteCarloBarostat(1*bar, simulation_temperature*kelvin))
    return sim_context(pdb, system, state, store_report_interval, nsteps, dcd_filename, platform, prop)

# System Config
platform, prop = set_compute_system("CUDA")
#platform, prop = set_compute_system("OpenCL")

# MD Protocol Control Flow
if restart:

    restart_number = currnt_restart_number +1
    #NVT
    if run_nvt_sim:
        print(f"Restarting simulation of {input_pdbfile}")
        pdb = PDBFile(f"nvt_md_{currnt_restart_number}.pdb")
        system, state = read_simulation_result(f"state_nvt_md_{currnt_restart_number}.xml", f"system_nvt_md_{currnt_restart_number}.xml")
        simulation = nvt_run(pdb, system, state, store_report_interval, equil_steps, f"nvt_md_{restart_number}", platform, prop)
        save_nvt_results = save_simulation_results(simulation, f"nvt_md_{restart_number}.pdb", f"state_nvt_md_{restart_number}.xml", f"system_nvt_md_{restart_number}.xml")
    if run_trex_sim:
        reporter = multistate.MultiStateReporter(f"REMD_{currnt_restart_number}.nc", open_mode="r")
        n_replicas = reporter.n_replicas
        n_iter_done = reporter.read_checkpoint_iterations().shape[0] -1
        n_iter_remaining = n_iter - n_iter_done
        sampler_states = reporter.read_sampler_states(n_iter_done)
        thermodynamic_states = reporter.read_thermodynamic_states()[0]
        ther_indices = [int(i) for i in reporter.read_replica_thermodynamic_states(iteration=n_iter_done)]
        sampler_states = [sampler_states[ther_indices.index(i)] for i in range(n_replicas)]
        move = mcmc.LangevinSplittingDynamicsMove(timestep=2*femtosecond, collision_rate=1/picosecond, n_steps=1000)
        simulation = multistate.ReplicaExchangeSampler(mcmc_moves=move, replica_mixing_scheme='swap-neighbors', number_of_iterations=n_iter_remaining)
        # Run the simulation and output the results in .nc file
        reporter = multistate.MultiStateReporter(f'REMD_{restart_number}.nc', checkpoint_interval=1)
        simulation.create(thermodynamic_states,sampler_states,storage=reporter)
        simulation.run()


# NEED TO RE-WRAP SIMULATION AFTER NVT AND T-REX 
else:
# Load parameters
    if run_nvt_sim:
        print("Loading parameters")
        params = read_params(path_to_topparstr)
        psf = read_top(path_to_psf)
        crd = read_crd(path_to_crdfile)
        psf = gen_box(psf, crd)
        system = psf.createSystem(params, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
        simulation = energy_minimization(psf, crd, system, platform, prop)
        simulation = rewrap(simulation)
        save_energy_min_results = save_simulation_results(simulation, "energy_min.pdb", "state_min.xml", "system_min.xml")
        system, state = read_simulation_result("state_min.xml","system_min.xml")
        pdb = PDBFile('energy_min.pdb')
        simulation = nvt_run(pdb, system, state, store_report_interval, equil_steps, f"nvt_md_{currnt_restart_number}", platform, prop)
        simulation = rewrap(simulation)
        save_nvt_equil_results = save_simulation_results(simulation, f"nvt_equil_{currnt_restart_number}.pdb", f"state_nvt_equil_{currnt_restart_number}.xml", f"system_nvt_equil_{currnt_restart_number}.xml")
        
    #T-REX MD
    if run_trex_sim:
        pdb = PDBFile(f'nvt_equil_{currnt_restart_number}.pdb')
        water_count = get_water_count(pdb)
        system, state = read_simulation_result( f"state_nvt_equil_{currnt_restart_number}.xml",f"system_nvt_equil_{currnt_restart_number}.xml")
        temperatures = temp_dict[water_count]* kelvin
        if min(temperatures) < 310* kelvin or max(temperatures) > 385* kelvin:
            print(f"min temperature found is {min(temperatures)} and max temperature found is {max(temperatures)}. These maybe beyond the range. Please check and resubmit.")
            sys.exit()
        n_replicas = len(temperatures)
        n_iterations = n_iter
        # Set the states, NPT ensemble
        thermodynamic_states = [ThermodynamicState(system=system, temperature=T, pressure=1*bar) for T in temperatures]
        sampler_state = SamplerState(positions=state.getPositions(), velocities=state.getVelocities(), box_vectors=state.getPeriodicBoxVectors())
        # Use Langevin dynamics integrator
        move = mcmc.LangevinSplittingDynamicsMove(timestep=2*femtosecond, collision_rate=1/picosecond, n_steps=1000)
        simulation = multistate.ReplicaExchangeSampler(mcmc_moves=move, replica_mixing_scheme='swap-neighbors', number_of_iterations=n_iterations)
        # Run the simulation and output the results in .nc file
        reporter = multistate.MultiStateReporter(f'REMD_{currnt_restart_number}.nc', checkpoint_interval=1)
        simulation.create(thermodynamic_states,sampler_state,storage=reporter)
        simulation.run()    
        #simulation = rewrap(simulation)
