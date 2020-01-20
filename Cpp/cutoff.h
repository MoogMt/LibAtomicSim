#ifndef CUTOFF_H
#define CUTOFF_H

#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>  
#include <stdlib.h>    
#include <math.h>

#include "utils.h"
#include "atom.h"
#include "cell.h"
#include "lut.h"

//================
// Cut-Off Matrix
//====================================================================
struct CutOffMatrix
{
  std::vector<double> matrix;
};
//====================================================================

//======================
// Read Cut Off Matrix
//=======================================================================================
bool readCutOff( std::ifstream & input , CutOffMatrix & com , std::vector<TypeLUT> & lut_list );
bool readCutOff( const std::string file , CutOffMatrix & com , std::vector<TypeLUT> & lut_list );
bool readCutOff(  std::ifstream & input , CutOffMatrix & com, AllTypeLUT & lut_list );
bool readCutOff( const std::string file , CutOffMatrix & com, AllTypeLUT & lut_list );
//=======================================================================================

//=============
// Get Cut Off 
//=======================================================================================
double getCutOff( const CutOffMatrix cut_off_matrix , const int i , const int j );
//=======================================================================================

#endif 
