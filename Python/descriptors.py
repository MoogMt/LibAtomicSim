#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 14:49:53 2019

@author: julienh with modifications from moogmt
"""

import tqdm
import numpy as np

from dscribe.descriptors import SOAP
from dscribe.descriptors import ACSF

# - SOAP
#=============================================================================#

default_sigma_SOAP  = 0.05  # Sigma for the gaussian density
default_cutoff_SOAP = 1.001 # Angstroms
default_nmax_SOAP   = 1     # nmax for the radial expansion
default_lmax_SOAP   = 0     # lmax for the angular expansion
default_sparse_SOAP = False # Sparse, avoid 0
        
def createDescriptorsSingleSOAP(data, species, sigma_SOAP, cutoff_SOAP, nmax_SOAP, lmax_SOAP, periodic,
                          sparse_SOAP=default_sparse_SOAP
                          ):

    # Initialize SOAP
    soap = SOAP( species=species, sigma=sigma_SOAP, periodic=periodic, rcut=cutoff_SOAP, nmax=nmax_SOAP, lmax=lmax_SOAP, sparse=sparse_SOAP )
    return soap.create( data )

def createDescriptorsAllSOAP( data, 
                              species, 
                              sigma_SOAP, 
                              cutoff_SOAP, 
                              nmax_SOAP, 
                              lmax_SOAP, 
                              periodic,
                              sparse_SOAP=default_sparse_SOAP
                            ):

    # Initialize SOAP
    soap = SOAP( species=species, sigma=sigma_SOAP, periodic=periodic, rcut=cutoff_SOAP, nmax=nmax_SOAP, lmax=lmax_SOAP, sparse=sparse_SOAP )
    
    # Compute number of features
    n_features = soap.get_number_of_features()
    n_atoms = np.shape(data[0])[0]
    n_steps = len(data)
    # Computing descriptors
    descriptors = np.empty( ( n_atoms, n_steps, n_features ),  dtype=object )
    for index_structure in tqdm.tqdm(range( n_steps )): 
        descriptors[:,index_structure,:] = soap.create( data[index_structure] )
    descriptors_ = []
    for atom in range( n_atoms ):
        descriptors_.append( descriptors[atom,:,:] )
    return descriptors_
#=============================================================================#



# ACSF
# TO BE REDONE AT A LATER POINT
#=============================================================================#
default_cutoff_ACSF = 0.                          # Cut-off for the ACSF function
default_n_acsf = 0                                # Number of ACSF functions per atoms
default_g2_params = np.zeros((default_n_acsf,2))  # Parameters for the g2 functions for the ACSF functions
default_g3_params = np.zeros((default_n_acsf,3))  # Parameters for the g3 functions for the ACSF functions
def createDescriptorsACSF(data, metadata,
                          cutoff_ACSF=default_cutoff_ACSF,
                          g2_params=default_g2_params,
                          g3_params=default_g3_params,
                          ):   
    # Updating metadata      
    
    # Prepping ACSF descriptor structure
    acsf = ACSF(species=metadata['species'],rcut=metadata['cutoff_acsf'],g2_params=metadata['g2_params'],g4_params=metadata['g3_params'])
    metadata['n_features'] = acsf.get_number_of_features()    
    # Computing descriptors
    descriptors = []  
    for index_atom in tqdm.tqdm(range(metadata['n_atoms'])):
        descriptors_loc = np.empty((np.shape(data)[0],metadata['n_features']))
        for index_structure in range(metadata['train_set_size']):
            descriptors_loc[index_structure,:] = acsf.create(data[index_structure],positions=[index_atom])
        descriptors.append( descriptors_loc )
    return metadata, descriptors
#=============================================================================#    

#from scipy.spatial.distance import cdist

#def create_data_NN(data,metadata):
#    local_dist_mat = np.empty((metadata['total_set_size'],metadata['n_atoms'],metadata['n_nearest_neigh']))
#    metadata['n_features'] = metadata['n_nearest_neigh']
#    for index in tqdm.tqdm(range(metadata['total_set_size'])):
#        dist_mat = cdist(data['structures'][index].get_positions(),data['structures'][index].get_positions())
#        for atom in range(metadata['n_atoms']):
#            if atom in [0,1]:
#                order = np.argsort(dist_mat[atom])[1:]
#                first_Hs = order[np.isin(order,range(2,7))][:metadata['n_nearest_neigh']]
#                first_Os = order[np.isin(order,range(0,2))][:metadata['n_nearest_neigh']]
#                local_dist_mat[index,atom] = dist_mat[atom,np.concatenate([first_Hs,first_Os])]
#            else:
#                order = np.argsort(dist_mat[atom])[1:]
#                first_Hs = order[np.isin(order,range(2,7))][:metadata['n_nearest_neigh']]
#                first_Os = order[np.isin(order,range(0,2))][:metadata['n_nearest_neigh']]
#                local_dist_mat[index,atom] = dist_mat[atom,np.concatenate([first_Hs,first_Os])]
#    return metadata, data.join(pd.DataFrame({'descriptor':list(local_dist_mat)}))
#    
    
    
    
    