#ifndef BONDS_H
#define BONDS_H

#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>  
#include <stdlib.h>    
#include <math.h>

#include "atom.h"
#include "contact_matrix.h"

//======
// BOND
//===============================
struct Bond
{
  int atom1_index, atom2_index;
  double strengh;
};
//===============================

//==========
// MAKE
//====================================================================================
Bond emptyBond( );
//-----------------------------------------------------------------------------------
Bond makeBond( int atom_index1 , Atom atom_index2 );
Bond makeBond( Atom atom1 , Atom atom2 );
//-----------------------------------------------------------------------------------
std::vector<Bond> makeBonds( std::vector<Atom> atom_list , double cut_off_radius );
std::vector<Bond> makeBonds( Contact_Matrix contact_matrix , double cut_off_radius );
//===================================================================================

#endif 
