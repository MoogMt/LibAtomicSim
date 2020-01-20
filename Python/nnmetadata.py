#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 12:15:45 2019

@author:  julienh with modification from moogmt

Contains functions for handling atomic simulation data
"""

import os
import ase
import numpy as np
import periodicTable as pT
import cpmd
import filexyz

from sklearn.preprocessing import StandardScaler  
from sklearn.decomposition import PCA


#==============================================================================
def handleOptionArg( input_label, default_value, metadata ):
    if not input_label in metadata :
        metadata[input_label] = default_value
    return metadata
def handleAllOptionaArg( input_label_list, default_values_list, metadata ):
    for i in range( len(input_label_list) ):
        metadata = handleOptionArg( input_label_list[i], default_values_list[i], metadata )
    return metadata
#==============================================================================

#==============================================================================
def getIndexExcluding( index_excluding, size_total ):
    return np.array(list(filter(lambda x : x not in index_excluding, np.arange(size_total))))
#==============================================================================

#============================================================================== 
def getDataFromIndex( input_, output_, indexs_chosen ):
    return input_[ indexs_chosen ], output_[ indexs_chosen ]
#------------------------------------------------------------------------------
def samplingExtremitiesFraction( fraction, energies, comp_time, nb_bins, comp_time_threshold ):
    if fraction > 1 :
        fraction = 1
    elif fraction < 0:
        fraction = 0
    if comp_time_threshold < 0:
        comp_time_threshold = 0
    nb_point=len(energies)
    # Chosen point will be 1, others 0
    choice_points=np.zeros( nb_point ,dtype=int)    
    # Keeping points that are above a given computational threshold
    choice_points[ np.nonzero( comp_time[ comp_time > comp_time_threshold ] )   ] = 1
    # Computing bins
    min_energy=energies.min()
    delta_energy=(energies.max()-min_energy)/nb_bins
    point_bins=np.array(np.round((energies-min_energy)/delta_energy,0),int)
    # Keeping the first and last bins 
    choice_points[ point_bins == 0 ] = 1
    choice_points[ point_bins == nb_bins-1 ] = 1 
    for point in range( nb_point ):
            if np.random.rand() > 1 - fraction:
                choice_points[point] += 1
    return np.nonzero(choice_points[ choice_points > 0 ])[0]
#------------------------------------------------------------------------------
def extractTrajectory( traj, chosen_index ):
    nb_train_size=len(chosen_index)
    structures = np.empty( nb_train_size, dtype=ase.atoms.Atoms )
    for i in range( nb_train_size ):
        structures[i] = traj[chosen_index[i]]
    return structures
#------------------------------------------------------------------------------
def choseDataRandomExclusion( traj, energies, test_size, index_exclude, replace ) :
    chosen_index=[]
    total_size=len(energies)
    if not replace:
        chosen_index = getIndexExcluding( index_exclude, total_size )
        chosen_index = chosen_index[ np.random.choice( chosen_index, size=test_size, replace=False ) ]
    else:
        chosen_index = np.random.choice( total_size , size=test_size, replace=True )
    energies=energies[chosen_index]
    structures = np.empty( test_size, dtype=ase.atoms.Atoms )
    for i in range( len(chosen_index) ) :
         structures[i] = traj[ chosen_index[i] ]
    return structures, energies
        
#==============================================================================

# Principal Component Analysis stuff
#==============================================================================
def makePCA( n_components ):
    return PCA(n_components=n_components)
#------------------------------------------------------------------------------
default_write_var = False
default_path_pcavar = "./"
def pcaVariance( input_, write_var=default_write_var,  path_pcavar=default_path_pcavar ):
    nb_point=len(input_[:])*len(input_[0])
    pca=PCA().fit( np.array(input_[:]).reshape(nb_point,np.shape(input_[0])[1]) )
    var=np.cumsum(pca.explained_variance_ratio_)
    if write_var :
        file_out=open(path_pcavar,"w")
        for i in range( np.shape(var)[0] ):
            file_out.write( str(i)+" "+str(var[i])+"\n" )
        file_out.close()
    return var
#------------------------------------------------------------------------------
default_verbose_pca = False
def pcaSelectBestParams( input_, n_pca, species, start_species, nb_element_species, n_features, verbose=default_verbose_pca ):
    pca=[]
    for specie in range(len(species) ):
        pca.append( PCA( n_components=n_pca ).
                   fit( input_[:,start_species[specie]:start_species[specie]+nb_element_species[specie],:].
                   reshape( input_[ :, start_species[specie]:start_species[specie]+nb_element_species[specie], : ].
                   shape[0]*nb_element_species[specie], n_features ) ) )
        if verbose: 
            print("Precision of new features ",specie," :",np.cumsum(pca[specie].explained_variance_ratio_)[-1],"\n" )        
    return pca
def pcaNFromVar( input_, fraction, species, start_species, nb_element_species, n_features, verbose=default_verbose_pca ):
    n_pca=np.zeros(2,dtype=int)
    for specie in range(len(species) ):
        pca=PCA( n_components=n_pca ).fit( input_[:,start_species[specie]:start_species[specie]+nb_element_species[specie],:].reshape( input_[ :, start_species[specie]:start_species[specie]+nb_element_species[specie], : ].shape[0]*nb_element_species[specie], n_features ) ) 
        variance=np.cumsum(pca[specie].explained_variance_ratio_)[-1]
        for i in range( len(variance) ):
             if variance[i] > fraction:
                 n_pca[specie] = i
                 break
    return n_pca
#==============================================================================

# Scalers - to be generalized
#==============================================================================
def createScaler( input_, species, start_species, nb_element_species, train_set_size, n_features ):
    # Create a scaler for the input using Scitkit Learn Standard Scaler
    scalers = []
    n_species = len(species)
    for specie in range( n_species ):
        scalers.append(StandardScaler())
        nb_atoms_total = nb_element_species[specie]*train_set_size
        start_specie = start_species[specie]
        end_specie   = start_specie + nb_element_species[specie]
        scalers[specie].fit(np.array(input_[start_specie:end_specie]).reshape(nb_atoms_total,n_features))    
    return scalers
#------------------------------------------------------------------------------
def applyScale( scalers, input_, species, start_species, nb_element_species, n_features ):
    n_species=len(species)
    # Applies a created scaler
    input_array=np.array(input_[:])
    for specie in range( n_species ):
        start_specie = start_species[specie]
        end_specie   = start_specie + nb_element_species[specie]
        nb_points = np.shape(input_array[:,:,:])[1]
        nb_atoms_total = nb_element_species[specie]*nb_points
        input_array[start_specie:end_specie,:,:] = scalers[specie].transform( input_array[start_specie:end_specie,:,:].reshape(nb_atoms_total, n_features ) ).reshape( nb_element_species[specie], nb_points, n_features )
    for i in range(np.shape(input_array)[0]):
        input_[i] = input_array[i,:,:]
    return input_
#------------------------------------------------------------------------------
#  Basic hand_made scalers
def scaleData( data ):
    min_data = data.min()
    range_data = data.max() - min_data
    return ( data - min_data)/range_data, min_data, range_data
#------------------------------------------------------------------------------
def deScaleData( data_scaled, min_data, range_data ):
    return data_scaled*range_data+min_data
#==============================================================================


# Handling Atoms stuff - may be moved at a later point
#==============================================================================
def getNbAtoms( atoms ):
    if type(atoms) == ase.atoms.Atoms :
        return len(atoms)
    elif type(atoms) == np.ndarray :
        return len(atoms)
    elif type(atoms) == list :
        if type(atoms[0]) == ase.atoms.Atoms :
            return len(atoms[0])
        elif type(atoms[0]) == np.ndarray :
            return len(atoms)
#------------------------------------------------------------------------------
def getSpecies( atoms ):
    if type(atoms) == list:
        atoms=atoms[0]
    types_=[]
    types_=np.append(types_,pT.z2Names(atoms.numbers[0]))
    for atom in range(getNbAtoms(atoms)):
        check = True
        for type_ in range(len(types_)):
            if pT.z2Names(atoms.numbers[atom]) == types_[type_]:
                check = False
        if check :
            types_=np.append(types_,pT.z2Names(atoms.numbers[atom]))
    return types_
#------------------------------------------------------------------------------
def getNbAtomsPerSpecies( atoms, species ):
    if type(atoms) == list:
        atoms=atoms[0]
    n_species = len(species)
    n_atoms  = len(atoms)
    nb_element_species = np.zeros( n_species, dtype=int )
    for specie in range( n_species ):
        for atom in range( n_atoms ):
            if species[specie] == pT.z2Names(atoms.numbers[atom]) :
                nb_element_species[specie] += 1
    return nb_element_species
#------------------------------------------------------------------------------
def sortAtomsUniq( atoms ):
    # Sorting by decreasing Z 
    # BEHOLD THE BUBLE SORT OF DEATH - Needs to be changed asap by something 
    # less ineficient - probably are python functions that can do stuff
    for i in range(len(atoms)):
        for j in range(i+1,len(atoms)):
            if atoms.numbers[i] < atoms.numbers[j]:
                store_z = atoms.numbers[i]
                store_positions = atoms.positions[i,:]
                atoms.numbers[i] = atoms.numbers[j]
                atoms.positions[i,:] = atoms.positions[j,:]
                atoms.numbers[j] = store_z
                atoms.positions[j,:] = store_positions
    return atoms
#------------------------------------------------------------------------------
def buildInput(data):
    input_=[]
    data_=np.stack(data["descriptor"].str[:].values)
    nb_atoms=np.shape(data_)[1]
    for atom in range(nb_atoms):
        input_.append(data_[:,atom,:])
    return input_
#------------------------------------------------------------------------------
def sortAtomsSpecie( atoms ) :
    if type(atoms) == list:
        for index in range(len(atoms)):
            atoms[index] = sortAtomsUniq(atoms[index])
        return atoms
    else:
        return sortAtomsUniq(atoms)
#------------------------------------------------------------------------------
def getStartSpecies( atoms, species ):
    if type(atoms) == list: 
        atoms=atoms[0]
    n_species=len(species)
    n_atoms=len(atoms)
    start_species = np.zeros( n_species, dtype=int )
    for specie in range( n_species ):
        for atom in range( n_atoms ):
            if species[specie] == pT.z2Names(atoms.numbers[atom]) :
                start_species[specie] = atom
                break
    return start_species
#==============================================================================

# Default physics params
default_n_atoms        = 1
default_n_species      = 1
default_species        = []
default_masses         = np.zeros(default_n_atoms)
default_pbc            = None
default_periodic       = False
default_total_size_set = 0
default_start_species  = np.zeros(default_n_atoms,dtype=int)
default_nb_element_species = np.zeros(default_n_atoms,dtype=int)
default_species_sorted = False
default_replace = False   # Whether or not we can select several times the same data point in the training set

# Default Descriptors parameters
default_descriptor   = 'None'

# Data treatment
default_PCA = False
default_pca_N = 0 
default_scale = False

# Default I/O settingss
default_path_to_import_model = ""        # If taking over from previous model, path where to save the NN
default_prefix = ""                      # prefix for the output files
default_suffix = ""                      # suffix for the output files

# Verbose
default_verbose = False
default_n_jobs=1


def buildMetaData( traj_file, energy_file, output_folder, 
                  n_atoms=default_n_atoms,
                  n_species=default_n_species,
                  species=default_species,
                  masses=default_masses,
                  pbc=default_pbc, 
                  periodic=default_periodic,
                  total_size_set = default_total_size_set,
                  start_species=default_start_species,
                  nb_element_species=default_nb_element_species,
                  species_sorted=default_species_sorted,
                  replace = default_replace,
                  descriptor=default_descriptor, 
 #                 neigh_lim=default_neigh_lim,
                  pca=default_PCA,
                  pca_n=default_pca_N,
                  scale=default_scale, 
                  path_to_import_model=default_path_to_import_model,
                  prefix=default_prefix,
                  suffix=default_suffix,
                  verbose=default_verbose,
                  n_jobs=default_n_jobs,
                  ):
    
    # Building Metadata Dic
    metadata={
            # Physics
            'n_atoms':n_atoms,
            'n_species':n_species,
            'species':species,
            'masses': masses,   
            'pbc': pbc,
            'periodic': periodic,
            'total_size_set': total_size_set,
            'start_species': start_species,
            'nb_element_species': nb_element_species,
            'species_sorted': species_sorted,
            'replace': replace,

             #Descriptor
            'descriptor_type': descriptor,    #Choice of descriptor
            'n_features': 1,
            
            # PCA or not PCA
            'pca_check': pca,             #False if no PCA, a number if PCA to select N first axis
            'pca_n': pca_n,
            'scaler': scale,             #False if no scaling, None if scaling
            # I/O        
            'traj_file': traj_file,
            'energy_file': energy_file,
            "path_2_model": path_to_import_model,    # To import from precedent NN, give datetime, or None
            'prefix_out': prefix,
            'suffix_out': suffix,
            "output_folder": output_folder,

            # Verbose
            'verbose': verbose,
            'n_jobs': n_jobs,

            }
    return metadata

def checkMetaDataIO( metadata, verbose ):
    if not os.path.isfile( metadata['traj_file'] ): 
        if verbose:
            print("Necessary input trajectory file does not exists! Exiting...")
        return False
    elif not os.path.isfile( metadata['energy_file'] ):
        if verbose:
            print("Necessary input Energy file does not exists! Exiting...")
        return False
    if not os.path.isdir( metadata['output_folder']):
        if verbose:
            print("Output folder does not exists, creating it...")
        os.mkdir(metadata['output_folder'])
    return True