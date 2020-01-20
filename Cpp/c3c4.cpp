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
  std::ofstream c2_angles_out("c2angles.dat");
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

  //-----------------
  // Reading Cut-Off
  //-------------------------------------------------------------------
  CutOffMatrix cut_off;
  if ( ! readCutOff( "cut_off.dat" , cut_off , lut_list ) )
    {
      return 1;
    }
  //-------------------------------------------------------------------

  //-------------------
  // Reading XYZ file
  //----------------------------------------------------
  while( readStepXYZfast( input , atom_list , lut_list, true, true ) )
    {
      if ( step % comp_step == 0 && step > start_step && step < end_step )
	{
	  // Makes the contact matrix
	  makeContactMatrix( cm_connection , cm_distance , atom_list, cell , cut_off , lut_list );
	  // Making molecules
	  std::vector<Molecule> molecules = makeMolecules( cm_connection );
	  // Calculating angles
	  for ( int i=0 ; i < molecules.size() ; i++  )
	    {
	      for ( int j=0 ; j < molecules[i].atom_index.size() ; j++ )
		{
		  std::vector<double> angles = getAngleAtom( cm_distance , molecules[i] , molecules[i].atom_index[j] );
		  if ( molecules[i].names[j] == "C" && angles.size() == 3)
		    {
		      std::vector<int> per_atoms_index = getBonded( molecules[i] , molecules[i].atom_index[j] );
		      int index_center = molecules[i].atom_index[j];
		      std::vector<double> position_center = getPosition( atom_list , index_center );
		      std::vector<double> position_atom1  = getPosition( atom_list , per_atoms_index[0] );
		      std::vector<double> position_atom2  = getPosition( atom_list , per_atoms_index[1] );
		      std::vector<double> position_atom3  = getPosition( atom_list , per_atoms_index[2] );
		      std::vector<double> vector1 = Difference( position_atom1 , position_atom2 );
		      std::vector<double> vector2 = Difference( position_atom1 , position_atom3 );
		      double dist = getDistanceFromPlan( vector1 , vector2 , position_center, position_atom1 );
		    }
		}
	    }
	  // Step making
	  std::cout << step << std::endl;
	}      
      step++;
     }
  //----------------------------------------------------
  
  //--------------------
  // Print other values
  //----------------------------------------------------
  std::cout << "Other C Angles: " << std::endl;
  for ( int i=0 ; i < c_others_nb.size() ; i++ )
    {
      std::cout << "value: " << c_others_nb[i] << std::endl;
    }
  //----------------------------------------------------
  std::cout << "Other O Angles: " << std::endl;
  for ( int i=0 ; i < o_others_nb.size() ; i++ )
    {
      std::cout << "value: " << o_others_nb[i] << std::endl;
    }
  //----------------------------------------------------
  
  //--------------
  //Closing fluxes
  //----------------------
  input.close();
  c2_angles_out.close();
  //----------------------
      
  return 0;
}
