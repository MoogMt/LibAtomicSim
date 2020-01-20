#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 22:03:35 2019

@author: moogmt
"""

import numpy as np

#========================
# Cell related functions
#========================================================
def minDistance( x_ , x0_, a_ ):
    dx_ = x_ - x0_;
    if dx_ > a_*0.5: 
        return dx_-a_
    elif dx_< - a_*0.5: 
        return dx_ + a_;
    else: 
        return dx_;
#========================================================
        
#========================================================
def addVector( r_, r0_ ):
    dr_ = np.zeros(( r_[:,0].size, 3 ));
    for i in range(r_[:,0].size):
        for j in range( 3 ):
            dr_[i,j] = r_[i,j] + r0_[i,j];
    return dr_;
#========================================================