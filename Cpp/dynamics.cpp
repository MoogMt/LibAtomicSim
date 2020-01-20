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
  //---------------------------------

  //--------
  // Output
  //--------------------------------------------------------------------------------
  std::ofstream C2_out ("C2_out.dat",  std::ios::out );
  std::ofstream C3_out ("C3_out.dat",  std::ios::out );
  std::ofstream C4_out ("C4_out.dat",  std::ios::out );
  std::ofstream O1_out ("O1_out.dat",  std::ios::out );
  std::ofstream O2_out ("O2_out.dat",  std::ios::out );
  std::ofstream changesC ("changesC.dat",  std::ios::out );
  std::ofstream changesO ("changesO.dat",  std::ios::out );
  std::ofstream changes ("changes_out.dat",  std::ios::out );
  //--------------------------------------------------------------------------------
  
  //----------------------
  // Physical parameters
  //--------------------------------------
  int step       = 1;  // Step counter
  int start_step = 5000; // Start step
  int end_step   = 1000000; // End step
  int comp_step  = 1; // Frequency of computation
  //--------------------------------------

  //---------------
  // Initializers
  //--------------------------------------------------
  AtomList  atom_list;  // Atoms in cell
  AllTypeLUT lut_list; // LUT for types
  ContactMatrix cm;    // Contact Matrix
  //--------------------------------------------------

  //--------------  
  // Atom Indexes
  //------------------------------
  std::vector<int> atom_indexesC;
  std::vector<int> atom_indexesO;
  //------------------------------
  
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

  //------------------------------
  std::vector<double> coordinanceC;
  std::vector<double> coordinanceO;
  std::vector<double> lifetime;
  std::vector<double> lifetime_C;
  std::vector<double> lifetime_C2;
  std::vector<double> lifetime_C3;
  std::vector<double> lifetime_C4;
  std::vector<double> lifetime_O;
  std::vector<double> lifetime_O1;
  std::vector<double> lifetime_O2;
  std::vector<double> changes_all;
  std::vector<double> changes_C;
  std::vector<double> changes_O;
  //------------------------------
  
  //-------------------
  // Reading XYZ file
  //----------------------------------------------------
  while( readStepXYZfast( input , atom_list , lut_list, true, true ) )
    {
      if ( step % comp_step == 0 && step >= start_step && step < end_step)
	{
	  cm  = makeContactMatrixSoft( atom_list , cell , cut_off , lut_list , 1.75 , 10 , 50 );
	  if ( step == start_step )
	    {
	      int nb_atoms = atom_list.names.size();
	      // Get index of species
	      atom_indexesC = getSpecieIndex( lut_list, "C" );
	      atom_indexesO = getSpecieIndex( lut_list, "O" );
	      // Initializing
	      lifetime.assign( nb_atoms , 0. );
	      coordinanceC.assign( atom_indexesC.size(), 0. );
	      coordinanceO.assign( atom_indexesO.size(), 0. );
	      // Compute initial stuff
	       for ( int i=0 ; i < atom_indexesC.size() ; i++ )
		{
		  coordinanceC[i] = cumSum( getAtomContact( cm , atom_indexesC[i] ) );
		}
	       for ( int i=0 ; i < atom_indexesO.size() ; i++ )
		{
		  coordinanceO[i] = cumSum( getAtomContact( cm , atom_indexesO[i] ) );
		}
	    }
	  else
	    {
	      int c_changes=0; int o_changes=0;
	      for ( int i=0 ; i < atom_indexesC.size() ; i++ )
		{
		  double sum = cumSum( getAtomContact( cm , atom_indexesC[i] ) );
		  if( ( sum - 2.5 )*( coordinanceC[atom_indexesC[i]] - 2.5 ) < 0  )
		    {
		      if ( sum - 2.5 < 0 )
			{
			  lifetime_C2.push_back( lifetime[i] );
			  lifetime[i] = 0;
			}
		      else
			{
			  lifetime_C3.push_back( lifetime[i] );
			  lifetime[i] = 0;
			}
		      c_changes++;
		    }
		  else if ( ( sum - 3.5 )*( coordinanceC[i] - 3.5 ) < 0  )
		    {
		      if ( sum - 3.5 < 0 )
			{
			  lifetime_C3.push_back( lifetime[i] );
			  lifetime[i] = 0;
			}
		      else
			{
			  lifetime_C4.push_back( lifetime[i] );
			  lifetime[i] = 0;
			}
		      c_changes++;
		    }
		  else
		    {
		      lifetime[atom_indexesC[i]]++;
		    }
		  coordinanceC[i] = sum;
		}
	      changesC << step << " " << c_changes << std::endl;
	      for ( int i=0 ; i < atom_indexesO.size() ; i++ )
		{
		  double sum = cumSum( getAtomContact( cm , atom_indexesO[i] ) );
		  if( ( sum - 1.5 )*( coordinanceO[atom_indexesO[i]] - 1.5 ) < 0  )
		    {
		      if ( sum - 1.5 < 0 )
			{
			  lifetime_O1.push_back( lifetime[i] );
			  lifetime[i] = 0;
			}
		      else
			{
			  lifetime_O2.push_back( lifetime[i] );
			  lifetime[i] = 0;
			}
		      o_changes++;
		    }
		  else
		    {
		      lifetime[atom_indexesO[i]]++;
		    }
		  coordinanceO[i] = sum;
		}
	      changesO << step << " " << o_changes << std::endl;
	      changes << step << " " << c_changes + o_changes << std::endl;
	    }
	  std::cout << step << std::endl;
	}
      step++;
    }
  //----------------------------------------------------

  int nb_box = 100;
  writeHistogram( C2_out , normalizeHistogram( makeRegularHistogram( lifetime_C2 , min(lifetime_C2) , max(lifetime_C2) , nb_box ) ) );
  if ( lifetime_C3.size() > 0 )
    {
      writeHistogram( C3_out , normalizeHistogram( makeRegularHistogram( lifetime_C3 , min(lifetime_C3) , max(lifetime_C3) , nb_box ) ) );
    }
  if ( lifetime_C4.size() > 0 )
    {
      writeHistogram( C4_out , normalizeHistogram( makeRegularHistogram( lifetime_C4 , min(lifetime_C4) , max(lifetime_C4) , nb_box ) ) );
    }
  if ( lifetime_O1.size() > 0 )
    {
      writeHistogram( O1_out , normalizeHistogram( makeRegularHistogram( lifetime_O1 , min(lifetime_O1) , max(lifetime_O1) , nb_box ) ) );
    }
  if ( lifetime_O2.size() > 0 )
    {
      writeHistogram( O2_out , normalizeHistogram( makeRegularHistogram( lifetime_O2 , min(lifetime_O2) , max(lifetime_O2) , nb_box ) ) );
    }
  
  // Close
  C2_out.close();
  C3_out.close();
  C4_out.close();
  O1_out.close();
  O2_out.close();
  changesC.close();
  changesO.close();
  changes.close();
  
  return 0;
}
