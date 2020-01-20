#ifndef XYZ_H
#define XYZ_H

#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <stdio.h>  
#include <stdlib.h>    
#include <math.h>

#include "atom.h"
#include "lut.h"

//=====================
// Get Number of steps
//========================================
int getStepXYZ( std::ifstream & file );
//========================================

//======
// READ
//==============================================================================================
bool readStepXYZfast( std::ifstream & file , AtomList & atoms , AllTypeLUT & lut_list , const bool same_type , const bool verbose );
bool readStepXYZ( std::ifstream & file , AtomList & atoms , AllTypeLUT & lut_list , const bool same_type = true , const bool verbose = true );
bool readStepXYZ( std::ifstream & file , AtomList & atoms , std::vector<TypeLUT> & lut_list , const bool same_type = true , const bool verbose = true );
bool readStepXYZ( std::ifstream & file , std::vector<Atom> & atoms , const std::vector<TypeLUT> & lut_list , const bool same_type = true , const bool verbose = true );
std::vector<Atom> readstepXYZ( std::ifstream & file );
//==============================================================================================

//=======
// WRITE
//===========================================================================================
void writeXYZ( std::ofstream & file , const std::vector<Atom> atom_list );
void writeXYZ( std::ofstream & file , const std::vector<Atom> atom_list , const int step );
void writeXYZ( std::ofstream & file , const AtomList atom_list );
void writeXYZ( std::ofstream & file , const AtomList atom_list , const int step );
//===========================================================================================
  
#endif 
