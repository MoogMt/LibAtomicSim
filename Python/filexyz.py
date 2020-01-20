#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 22:03:09 2019

@author: moogmt

Functions that deals with XYZ file 
Mostly built-in function in the ASE library 
# Remember to add proper ASE citations...
"""

import ase.io
import numpy as np

def getNbStep( file_path ):
    return len(ase.io.read(file_path),index=':')

def read( file_path ):
    return ase.io.read(file_path,index=':')

def readPbcCubic( file_path, a ):
    traj=read( file_path )
    for i in range(len(traj)):
        traj[i].pbc=True
        traj[i].set_cell([a, a, a])
    return traj

def readAsArray( file_path ):
    traj=ase.io.read(file_path,index=':')
    traj_ = np.zeros( len(traj),3)
    for i in range( len(traj)):
        traj_[i,:] = traj[0].positions
    return 
