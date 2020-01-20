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
#include "lut.h"
#include "molecules.h"
//-------------------------

//================
// MAIN PROGRAM
//=====================================================================
int main( void )
{
  //--------
  // Input
  //---------------------------------
  std::ifstream input("TRAJEC.xyz");
  //--------------------------------

  //--------
  // Output
  //-------------------------------------
  std::ofstream mols_out("mols_out.dat");
   //-------------------------------------
  
  //----------------------
  // Physical parameters
  //--------------------------------------
  int step       = 1;  // Step counter
  int start_step = 2000; // Start step
  int end_step   = 120000;
  int comp_step  = 1; // Frequency of computation
  //--------------------------------------

  //---------------
  // Initializers
  //--------------------------------------------------
  AtomList  atom_list;  // Atoms in cell
  AllTypeLUT lut_list; // LUT for types
  ContactMatrix cm_connection;    // Contact Matrix
  ContactMatrix cm_distance;    // Contact Matrix
  //--------------------------------------------------

  //--------------------
  // Reading Cell File
  //-------------------------------------------------------------------
  Cell cell;
  if ( ! readParamCell( "cell.param" , cell ) )
    {
      return 1;
    }
  //-------------------------------------------------------------------


  std::vector<MolSig> mol_sig;
  std::vector<int> nb_sig;

  //-----------------
  // Reading Cut-Off
  //-------------------------------------------------------------------
  CutOffMatrix cut_off;
  if ( ! readCutOff( "cut_off.dat" , cut_off , lut_list ) )
    {
      return 1;
    }
  //-------------------------------------------------------------------

  // Reading XYZ file
  //----------------------------------------------------
  while( readStepXYZfast( input , atom_list , lut_list, true, true ) )
    {
      if ( step % comp_step == 0 && step > start_step && step < end_step )
	{
	  //--------------------------
	  // Makes the contact matrix
	  //----------------------------------------------------------------
	  makeContactMatrix( cm_connection , cm_distance , atom_list, cell , cut_off , lut_list );
	  //-----------------------------------------------------------------

	  //------------------
	  // Making molecules
	  //-----------------------------------------------------------------
	  std::vector<Molecule> molecules = makeMolecules( cm_connection );
	  //-----------------------------------------------------------------

	  //-----------------------------
	  // Computing and updating Sig
	  //----------------------------------------
	  updateSigs( molecules, mol_sig , nb_sig );
	  //----------------------------------------

	  //-----
	  // Step
	  //------------------------------
	  std::cout << step << std::endl;
	  //------------------------------
	}      
      step++;
     }
  //----------------------------------------------------
  
  //--------------
  //Closing fluxes
  //----------------------
  input.close();
  mols_out.close();
  //----------------------
      
  return 0;
}
