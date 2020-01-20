#ifndef MOLECULES_H
#define MOLECULES_H

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
#include "cutoff.h"
#include "contact_matrix.h"
#include "lut.h"
#include "bonds.h"

//==========
// MOLECULE
//=============================================
struct MoleculeBasic
{
  std::vector<std::string> names;
  std::vector<int> atom_index;
};
//---------------------------------------
struct Molecule
{
  std::vector<std::string> names;
  std::vector<int> atom_index;
  std::vector<Bond> bonds;
};
//=============================================

//======
// MAKE
//=============================================================================================
MoleculeBasic startMoleculeBasic( const std::string name_atom , const int index_atom );
std::vector<MoleculeBasic> makeMoleculesBasic( const ContactMatrix & cm );
//----------------------------------------------------------------------------------
Molecule startMolecule( const std::string name_atom , const int index_atom );
std::vector<Molecule> makeMolecules( const ContactMatrix & cm );
//============================================================================================

//=======
// PRINT
//============================================================================
void printMolecule( const MoleculeBasic molecule );
void printMoleculeSize( const MoleculeBasic molecule, const bool toline );
void printMolecules( const std::vector<MoleculeBasic> molecules );
void printMoleculesSize( const std::vector<MoleculeBasic> molecules );
//============================================================================

//========
// Bonded
//=================================================================================
std::vector<int> getBonded( Molecule molecule , int atom_index );
//=================================================================================

//========
// Angles 
//============================================================================
std::vector<double> getAngleAtom( ContactMatrix cm , Molecule molecule , int atom_index );
//============================================================================

//=========
// SIGNALS
//===============================================================================================
//int sigExists( MolSig signal , std::vector<MolSig> mol_signals )
//void UpdateSigs( std::vector<Molecule> molecules , std::vector<MolSig> & mol_signals , std::vector<int> & occurences )
//===============================================================================================
#endif
