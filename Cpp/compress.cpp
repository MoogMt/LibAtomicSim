//---------------
// External Libs
//-------------------
#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>  
#include <stdlib.h>    
#include <math.h>
//---------------------

//----------------
// Internal Files
//-------------------------
#include "utils.h"
#include "atom.h"
#include "cell.h"
#include "contact_matrix.h"
#include "xyz.h"
#include "histogram.h"
#include "pdb.h"
#include "sim.h"
//-------------------------

//================
// MAIN PROGRAM
//=====================================================================
int main(void)
{
  //--------
  // Input
  //---------------------------------
  std::ifstream input("input.xyz");
  //---------------------------------

  //--------
  // Output
  //-------------------------------------------------------------
  std::ofstream C_position( "C_position.cpmd" , std::ios::out );
  std::ofstream O_position( "O_position.cpmd" , std::ios::out );
  std::ofstream xyz( "position.xyz" , std::ios::out );
  //-------------------------------------------------------------

  //----
  // As
  //-------------------------------------------------------------
  double a_prev;
  std::cout << "Give the present a: ";
  std::cin >> a_prev;
  double b_prev;
  std::cout << "Give the present b: ";
  std::cin >> b_prev;
  double c_prev;
  std::cout << "Give the present c: ";
  std::cin >> c_prev;
  double a_new;
  std::cout << "Give the new a: ";
  std::cin >> a_new;
  double b_new;
  std::cout << "Give the new b: ";
  std::cin >> b_new;
  double c_new;
  std::cout << "Give the new c: ";
  std::cin >> c_new;
  //-------------------------------------------------------------
  
  //----------------------
  // Physical parameters
  //----------------------------------------------------------------
  double frac_a = a_new/a_prev;
  double frac_b = b_new/b_prev;
  double frac_c = c_new/c_prev;
  //----------------------------------------------------------------

  //---------------
  // Initializers
  //----------------------------------------------
  AtomList atom_list;
  AllTypeLUT lut_list;
  //----------------------------------------------

  //---------------
  // Reading cell
  //---------------------------------
  Cell cell = { a_prev , b_prev , c_prev , 90 , 90 , 90 };
  //---------------------------------

  //-------------------
  // Reading XYZ file
  //----------------------------------------------------
  if ( readStepXYZfast( input , atom_list , lut_list , false , false) )
    {
      wrapPBC( atom_list , cell) ;
      compressCell( atom_list , cell , frac_a , frac_b , frac_c );
      writeXYZ( xyz , atom_list );
      writePositions( C_position , atom_list , "C" );
      writePositions( O_position , atom_list , "O" );
    }
  //----------------------------------------------------

  //Closing fluxes
  //----------------------
  input.close();
  C_position.close();
  O_position.close();
  xyz.close();
  //----------------------
  
  return 0;
}
