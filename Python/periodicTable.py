#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 10:58:58 2020

@author: moogmt
"""

zNames=[ "H",                                                                                                                                                                                     "He",
        "Li", "Be",                                                                                                                                                  "B",  "C",  "N",  "O",  "F", "Ne", 
        "Na", "Mg",                                                                                                                                                 "Al", "Si",  "P",  "S", "Cl", "Ar",
         "K", "Ca",                                                                                     "Sc", "Ti",  "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", 
        "Rb", "Sr",                                                                                      "Y", "Zr", "Nb", "Mo", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sb", "Te",  "I", "Xe",
        "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta",  "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Ti", "Pb", "Bi", "Po", "At", "Rn",
        "Fr", "Ra", "Ac", "Th", "Pa",  "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"]

def z2Names( z ):
    return zNames[z-1]
        
def names2Z( name ):
    for i in range(len(zNames)):
        if name == zNames[i]:
            return i+1
        