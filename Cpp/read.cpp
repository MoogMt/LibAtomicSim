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
  std::ofstream c3_angles_out("c3angles.dat");
  std::ofstream c4_angles_out("c4angles.dat");
  std::ofstream o2_angles_out("o2angles.dat");
  std::ofstream o3_angles_out("o3angles.dat");
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

  //------------
  // Histograms
  //---------------------------------------------------------------
  // Technical values
  double hist_start = 0;
  double hist_end   = 180;
  int nb_box        = 1800;
  //---------------------------------------------------------------
  std::vector<Bin> c2_angles_hist; std::vector<double> c2_angles;
  std::vector<Bin> c3_angles_hist; std::vector<double> c3_angles;
  std::vector<Bin> c4_angles_hist; std::vector<double> c4_angles;
  std::vector<Bin> o2_angles_hist; std::vector<double> o2_angles;
  std::vector<Bin> o3_angles_hist; std::vector<double> o3_angles;
  std::vector<int> c_others_nb;    std::vector<int>  o_others_nb;
  //---------------------------------------------------------------
  
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
		  if ( angles.size() == 0 ) continue;
		  if ( molecules[i].names[j] == "C" )
		    {
		      if ( angles.size() == 1 )   appendVector( c2_angles , angles );
		      else if ( angles.size() == 3 )  appendVector( c3_angles , angles ); 
		      else if ( angles.size() == 6 )  appendVector( c4_angles , angles );
		      else c_others_nb.push_back( angles.size() );
		    }
		  else
		    {
		      if ( angles.size() == 1 ) appendVector( o2_angles , angles );
		      else if ( angles.size() == 3 )  appendVector( o3_angles , angles ); 
		      else o_others_nb.push_back( angles.size() );
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
  // Making histograms
  //----------------------------------------------------
  c2_angles_hist = makeRegularHistogram( c2_angles , hist_start , hist_end , nb_box );
  writeHistogram( c2_angles_out , normalizeHistogram( c2_angles_hist ) );
  c3_angles_hist = makeRegularHistogram( c3_angles , hist_start , hist_end , nb_box );
  writeHistogram( c3_angles_out , normalizeHistogram( c3_angles_hist ) );
  c4_angles_hist = makeRegularHistogram( c4_angles , hist_start , hist_end , nb_box );
  writeHistogram( c4_angles_out , normalizeHistogram( c4_angles_hist ) );
  o2_angles_hist = makeRegularHistogram( o2_angles , hist_start , hist_end , nb_box );
  writeHistogram( o2_angles_out , normalizeHistogram( o2_angles_hist ) );
  o3_angles_hist = makeRegularHistogram( o3_angles , hist_start , hist_end , nb_box );
  writeHistogram( o3_angles_out , normalizeHistogram( o3_angles_hist ) );
  //----------------------------------------------------

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
  c3_angles_out.close();
  c4_angles_out.close();
  o2_angles_out.close();
  o3_angles_out.close();
  //----------------------
      
  return 0;
}
