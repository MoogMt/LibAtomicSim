<#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 16 18:06:47 2020

@author: mathieumoog
"""

import cpmd
import os
import numpy as np
from dscribe.descriptors import SOAP
from ase import Atoms
from PIL import Image

def array2png( data, file_out ):
    new_im = Image.fromarray( data )
    new_im = new_im.convert( 'RGB' )
    new_im.save( file_out )
    return True

# Folder set up
co2_path_base="/Users/mathieumoog/Documents/CO2/"
folder_img = "/Users/mathieumoog/Documents/ImgCO2/"

# Sim parameters
volume      = 9.8
temperature = 2000 
run_nb      = 1

# Target folder
target_path  = co2_path_base 
target_path += str(volume) + "/" 
target_path += str(temperature) + "K/" 
target_path += str(run_nb) + "-run/"


if not os.path.exists( target_path ):
    print("Find does not exists, processing.")


# SOAP Setting    
species = [6,8]
sigma_soap = 2.5
periodic_soap = True
cutoff_soap = 5.0
n_soap = 4
l_soap = 4
sparse_soap = False
rbf_soap = 'gto'
soap = SOAP( 
            species  = species, 
            sigma    = sigma_soap, 
            periodic = periodic_soap, 
            rcut     = cutoff_soap, 
            nmax     = n_soap, 
            lmax     = l_soap,
            rbf      = rbf_soap,
           )
size_n = int(n_soap*(n_soap+1)/2)
size_l = l_soap + 1 

# Structure parameters
nb_atoms = 96
nbC = 32
nbO = 64
c_index_base = np.arange(0,nbC)
o_index_base = np.arange(nbC,nbC+nbO)
n_species=2
# Number of SOAP blocks
n_blocks = int(n_species*( n_species + 1 )/2)

# RGB parameters
max_pix=255
min_pix=0

# Reading Trajectory
correctness = True # Use only sane FTRAJECTORY files
ftraj_path = target_path + "FTRAJECTORY_db"
positions, velocities, forces = cpmd.readFtraj( ftraj_path, correctness )
velocities = [] # We don't use velocities and forces
forces     = []
nb_steps = np.shape( positions )[0] # nb of simulation step
# Conversion from Bohr to Angstrom 
bohr2Ang = 0.529177
positions=positions*bohr2Ang
# Wrap
positions %= volume

# Index of atoms
atoms_species = np.zeros( nb_atoms )
for carbons in c_index_base:
    atoms_species[carbons] = 6
for oxygen in o_index_base:
    atoms_species[oxygen] = 8

# Image Creation parameter ( used for the carbon, to get a square image)
nx=8
ny=4

# Creating ASE Atoms Object
traj = np.zeros((nb_steps),dtype=Atoms)
c_indexs = []
for step in range(nb_steps):
    traj[step] = Atoms( numbers=atoms_species, 
                        positions=positions[step,:,:],
                        cell=[volume,volume,volume], 
                        pbc=[ True, True, True ] 
                        )
    c_indexs.append( c_index_base )
    
# At that point positions are loaded in traj
positions=[]

# Compute whole soap traj
C_soap = soap.create( traj, positions=c_indexs )

# Remove traj to preserve memory
traj=[]

# Rescale parameter
rescale_specie = False

# Loop over steps
for step in range( nb_steps ):
    # Separating blocks of the SOAP vector
    data_raw = np.zeros(( size_n*4, size_l*8, 3 ))
    for i in range(ny):
        for j in range(nx):
            for block in range(n_blocks):
                start_block = block*size_n*size_l 
                end_block   = ( block + 1 )*size_n*size_l 
                data_raw[ i*size_n:(i+1)*size_n, j*size_l:(j+1)*size_l, block ] = np.reshape( C_soap[ (i*nx + j) + nbC*step , start_block:end_block ], ( size_n, size_l ) )

    # Min Max Scaling
    if rescale_specie: 
        for block in range(n_blocks):
            data_raw[:,:,block] = ( data_raw[:,:,block] - np.min(data_raw[:,:,block]) )/( np.max(data_raw[:,:,block])-np.min(data_raw[:,:,block]) )
    else:
        data_raw = ( data_raw - np.min(data_raw) )/( np.max(data_raw)-np.min(data_raw) )
        
    # Creating RGB Data
    data_pic = np.uint8(data_raw*( max_pix - min_pix ) + min_pix)

    step_print=str(step)
    folder_CC = folder_img + "CC/"
    _ = array2png( data_pic[:,:,0], str( folder_CC + step_print + ".png" ) )
    folder_CO = folder_img + "CO/"
    _ = array2png( data_pic[:,:,1], str( folder_CO + step_print + ".png" ) )
    folder_OO = folder_img + "OO/"
    _ = array2png( data_pic[:,:,2], str( folder_OO + step_print + ".png" ) )
    folder_all = folder_img + "All/"
    _ = array2png( data_pic[:,:,:], str( folder_all + step_print + ".png" ) )

# To make into a movie:
#ffmpeg -start_number 1 -i '%1d.png' -vcodec libx264 -crf 22 output.mp4
    
    
    
    
    
