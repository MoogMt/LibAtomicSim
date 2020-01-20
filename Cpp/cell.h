#ifndef CELL_H
#define CELL_H

#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>  
#include <stdlib.h>    
#include <math.h>
#include <sstream>
#include "utils.h"
#include "atom.h"

//======
// CELL
//===========================
struct Cell
{
  double a, b, c;             // lengths parameters
  double alpha, beta, gamma;  // angular parameters
};
//============================

//=====
// PBC
//==========================================================
double backIn( double x, double a );
std::vector<double> backIn ( std::vector<double> x , double a );
std::vector<double> wrapPBC( std::vector<double> vector_input, Cell cell );
Atom wrapPBC(Atom atom_in, Cell box );
std::vector<Atom> wrapPBC( std::vector<Atom> atoms , Cell cell );
void wrapPBC( AtomList & atom_list , int i , Cell cell );
void wrapPBC( AtomList & atom_list , Cell cell );
std::vector<Atom> pbcImages( Atom atom, Cell box );
void unwrap( AtomList & atom_list , AtomList & atom_list0 );
std::vector<double> getMinImage( std::vector<double> x , std::vector<double> x_ref , double cell_length );
std::vector<double> getMinImage( AtomList atom_list , Cell cell , int atom_center , int atom_target );
//==========================================================

//=======================
// Calculating Distances
//===============================================================================================
double distAtoms1D( double x1, double x2, double a );
double distanceAtoms( std::vector<Atom> atoms, int i, int j, Cell box , bool wrap = true , bool sqrt_test = true );
double distanceAtomsSq( AtomList & atom_list ,  int i , int  j , Cell cell , bool wrap=true );
//===============================================================================================

//===========
// Change Box
//============================================================================================
void compressCell( Cell & cell , double frac_a , double frac_b , double frac_c );
void compressAtoms( AtomList & atom_list, double frac_a , double frac_b , double frac_c );
void compressCell( AtomList & atom_list , Cell & cell , double frac_a , double frac_b , double c );
//=============================================================================================

//====
// IO
//=============================================================================================
// Writting Atom distances
void writeAtomDistances( std::ofstream & file , std::vector<Atom> atom_list , std::vector<int> atom_index, Cell box);
// Reading Cell Param File
bool readParamCellStep( std::ofstream & file , Cell & cell );
bool readParamCell( std::string file_name , Cell & cell );
//=============================================================================================

//==========
// Pressure
//=========================================================================
bool readPressureCell( std::ifstream & input , double & pressure );
bool readPressure( std::ifstream & input , double & pressure );
//=========================================================================
  
#endif 
