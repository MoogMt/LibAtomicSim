#ifndef SIM_H
#define SIM_H

#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>  
#include <stdlib.h>    
#include <math.h>

#include "atom.h"
#include "cell.h"
#include "utils.h"
#include "lut.h"
#include "xyz.h"

//=====
// SIM
//===========================
struct Sim
{
  std::vector<Atom> atoms;
  Cell cell;
};
//===========================

//==========
// Function
//===========================
Sim emptySim();
Sim compressBox( Sim sim_set , double frac_a , double frac_b , double frac_c );
//================================================================================================

//=======================
// Diffusion Coefficient
//================================================================================================
void computeDiff( std::ofstream & output , std::ifstream & input , int comp_step , int start_step , int end_step , AtomList atom_list , AllTypeLUT lut_list , Cell cell );
//================================================================================================

#endif 
