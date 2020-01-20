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