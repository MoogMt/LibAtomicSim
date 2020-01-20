#ifndef LUT_H
#define LUT_H

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

//=====
// LUT
//=============================
struct TypeLUT
{
  std::string name;
  int index;
  std::vector<int> atom_index;
};
//==============================
struct AllTypeLUT
{
  std::vector<int> type_index;
  std::vector<std::string> type_name;
  std::vector<TypeLUT> types;
};
//===============================

//=======
// Tests
//=============================================================================================
bool testType ( const AllTypeLUT lut_all, const int atom_index ,  const std::string specie );
bool testType ( const AllTypeLUT lut_all, const int atom_index ,  const int type_index );
bool testType ( const TypeLUT lut , const std::string specie );
bool testType ( const TypeLUT lut , const int index );
//--------------------------------------------------------------------------
bool typeExists( const std::vector<TypeLUT> list , const std::string specie );
bool typeExists( const std::vector<TypeLUT> list , const int index );
//=============================================================================================

//=============
// Check Names
//=============================================================================================
bool checkNames( std::vector<std::string> names );
//=============================================================================================

//========
// Modify
//=============================================================================================
void addAtom2LUT( std::vector<TypeLUT>  & list , const std::string name , const int index  );
void addAtom2LUT( std::vector<TypeLUT> & list , const Atom atom );
void addAtom2LUT( AllTypeLUT & list , const std::string name , const int index  );
void addAtom2LUT( AllTypeLUT & list , const Atom atom );
//=============================================================================================

//========
// Create
//=============================================================================================
TypeLUT makeLUT( const std::string name , const int index );
TypeLUT makeLUT( const std::string name , const int index , const std::vector<int> atom_index );
void makeLUT( std::vector<TypeLUT> & lut_list , const std::vector<std::string> names );
std::vector<TypeLUT> makeLUT( const std::vector<Atom> atoms );
void makeLUT( AllTypeLUT & lut_list , const std::vector<std::string> names );
//=============================================================================================

//============
// Get Indexes
//=============================================================================================
std::vector<int> getSpecieIndex( const AllTypeLUT & lut_maj, std::string specie );
//=============================================================================================

#endif
