#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 16 18:06:47 2020

@author: mathieumoog
"""

import numpy as np
from dscribe.descriptors import SOAP
from ase import Atoms
import matplotlib.pyplot as plt

from quippy import descriptors

# SOAP Setting    
species = [6,8]
sigma_soap = 0.05
periodic_soap = True
cutoff_soap = 6
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


volume=10

atoms_species=[6,8,8,6,8,8]

dist = 1.2
dist2 = 2

positions=np.ones((6,3),dtype=float)*volume*0.5
#Carbon is at the center of the cell
positions[1,0] += dist
positions[2,0] -= dist
positions[3:6,:] = positions[0:3,:] 
positions[3:6,1] += dist2

test = Atoms( numbers=atoms_species, 
                        positions=positions[:,:],
                        cell=[volume,volume,volume], 
                        pbc=[ True, True, True ] 
                        )


desc = descriptors.Descriptor("soap cutoff=3.5 l_max=4 n_max=4 atom_sigma=0.05 n_Z=2 Z={6 8} ")

data_desc = desc.calc(test)

data = data_desc["data"]

plt.plot(data[0,:])
plt.plot(data[1,:])
plt.plot(data[2,:])
plt.plot(data[3,:])
plt.plot(data[4,:])
plt.plot(data[5,:])

data.shape

n_species=2
n_block = int(n_species*( n_species + 1 )/2)

C_soap = soap.create( test)

data = np.zeros(( size_l,size_n, n_block ), dtype=float)

for block in range(n_block):
    start_block = block*size_n*size_l 
    end_block   = ( block + 1 )*size_n*size_l 
    data[:,:,block] = np.reshape( C_soap[ 0, start_block:end_block ], ( size_l, size_n ) )

plt.imshow(data[:,:,0] ) 
plt.show()

plt.imshow(data[:,:,1] )
plt.show()

plt.imshow(data[:,:,2] )
plt.show()













