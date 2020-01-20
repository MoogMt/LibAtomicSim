#ifndef ATOMS_H
#define ATOMS_H

#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>  
#include <stdlib.h>    
#include <math.h>

#include "utils.h"

//======
// ATOM
//==============================
struct Atom
{
  std::string name;    // name
  double x, y, z; // position
  int index;      // index
};
//==============================

//===========
// ATOM LIST
//================================
struct AtomList
{
  std::vector<std::string> names;
  std::vector<double> x, y, z;
  std::vector<int> index;
};
//================================

//=========
// Position
//=============================================================================================
std::vector<double> getPosition( AtomList & atoms , int index );
//=============================================================================================

//===========
// Distances
//=============================================================================================
double distanceAtoms(Atom i, Atom j);
std::vector<double> distanceFromPoint( AtomList atom_list , std::vector<double> point );
//=============================================================================================

//=======
// MOVE 
//=============================================================================================
std::vector<Atom> compressAtoms( std::vector<Atom> atoms, double frac_a , double frac_b , double frac_c );
//=============================================================================================

//=====
// IO
//=============================================================================================
// WRITE
//-----------------------------------------------------------------------------------------
void writePositions( std::ofstream & file , std::vector<Atom> atoms, std::string specie );
void writePositions( std::ofstream & file , AtomList atom_list , std::string specie );
//=============================================================================================

#endif 
