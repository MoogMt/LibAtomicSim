#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 27 15:34:48 2019

@author: moogmt
"""

import numpy as np

# Column code for ENERGIES file
cpmd_temperature_col=2
cpmd_pot_energy_col=3
cpmd_tot_energy_col=4
cpmd_msd_col=6
cpmd_scf_comptime=7

def getNbLineEnergies(file_path):
    f=open(file_path,"r") 
    lines=f.readlines() 
    f.close()
    return len(lines)

def readPotEnergy(file_path):
    nb_point=getNbLineEnergies(file_path)
    energies=np.zeros(nb_point)
    f=open(file_path,"r")
    for i in range (nb_point):
        energies[i] = f.readline().split()[cpmd_pot_energy_col] # Read Kohn-Sham energies (in Ry)
    f.close()
    return energies

def readEnergiesFile(file_path):
    nb_point=getNbLineEnergies(file_path)
    data=np.zeros((nb_point,7))
    f=open(file_path,"r")
    for i in range (nb_point):
        data[i,:] = f.readline().split()[1:] # Read all data except first column (time)
    f.close()
    return data

def extractTemperature(data):
    return data[:,cpmd_tot_energy_col-1]

def extractPotentialEnergy(data):
    return data[:,cpmd_pot_energy_col-1]

def extractTotalEnergy(data):
    return data[:,cpmd_tot_energy_col-1]

def extractMSD(data):
    return data[:,cpmd_msd_col-1]

def extractSCFcomputationTime(data):
    return data[:,cpmd_scf_comptime-1]


def getNbStepAtomLinesFtraj( file_input ):
    nb_atoms=0
    nb_lines=0
    with open( file_input ) as handle_in:
        for line in handle_in:
            keys = line.rsplit()
            if keys[0] != "<<<<<<" :
                if keys[0] == "1" :
                    nb_atoms+=1
                nb_lines+=1
    return int(nb_lines/nb_atoms), nb_atoms, nb_lines

def readFtraj( file_input, correctness ):
    nb_step, nb_atoms, nb_lines = getNbStepAtomLinesFtraj( file_input )
    if nb_lines % nb_atoms != 0 and not correctness:
        return False
    # Init vectors
    positions  = np.zeros( ( nb_step, nb_atoms, 3 ), dtype=float )
    velocities = np.zeros( ( nb_step, nb_atoms, 3 ), dtype=float )
    forces     = np.zeros( ( nb_step, nb_atoms, 3 ), dtype=float )        
    # Open file
    handle_in = open( file_input )
    # Loop over steps
    for step in range(nb_step):
        # Loop over atoms
        for atom in range(nb_atoms):
            # Read line
            all_line = handle_in.readline()
            # Parse Line
            line = all_line.rsplit()
            # Remove lines with issues
            if line[0] != "<<<<<<" :
                # Loop over dimensions
                for i in range(3):
                    # Positions
                    positions[  step, atom, i ] = line[ i + 1 ]
                    # Velocities
                    velocities[ step, atom, i ] = line[ i + 4 ]
                    # Forces
                    forces[     step, atom, i ] = line[ i + 7 ]
    handle_in.close()
    return positions, velocities, forces



