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
  std::ofstream c2_ratio_out("c2ratio.dat");
  std::ofstream c3_ratio_out("c3ratio.dat");
  std::ofstream c4_ratio_out("c4ratio.dat");
  std::ofstream c2_out("c2.dat");
  std::ofstream c3_out("c3.dat");
  std::ofstream c4_out("c4.dat");
   //-------------------------------------
  
  //----------------------
  // Physical parameters
  //--------------------------------------
  int step       = 1;  // Step counter
  int start_step = 1; // Start step
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
  std::vector<int> c_others_nb;
  std::vector<Bin> c2_ratio_hist ; std::vector<double> c2_ratio;
  std::vector<Bin> c3_ratio_hist ; std::vector<double> c3_ratio;
  std::vector<Bin> c4_ratio_hist ; std::vector<double> c4_ratio;
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
	  // Loop over all molecules
	  for ( int i=0 ; i < molecules.size() ; i++  )
	    {
	      // Loop over all atoms in a molecule
	      for ( int j=0 ; j < molecules[i].atom_index.size() ; j++ )
		{
		  // Computing all angles around an atom
		  std::vector<double> angles = getAngleAtom( cm_distance , molecules[i] , molecules[i].atom_index[j] );
		  // If there is no angles, we move on...
		  if ( angles.size() == 0 ) continue;
		  // We care only about C...
		  if ( molecules[i].names[j] == "C" )
		    {
		      if ( angles.size() == 1 )   appendVector( c2_angles , angles ); 
		      else if ( angles.size() == 3 )
			{
			  appendVector( c3_angles , angles );
			}
		      else if ( angles.size() == 6 )  appendVector( c4_angles , angles );  
		      else
			{
			  std::cout << "nb_angle: " << angles.size() << " ";
			  std::cout << "C_index: " <<  molecules[i].atom_index[j] << std::endl;
			  c_others_nb.push_back( angles.size() );
			}
		    }
		}
	    }
	  // Computing ratios
	  double angles_total = c2_angles.size() + c3_angles.size() + c4_angles.size();
	  //---------------
	  c2_out << step << " " << (double)(c2_angles.size())/angles_total << std::endl;
	  c3_out << step << " " << (double)(c3_angles.size())/angles_total << std::endl;
	  c4_out << step << " " << (double)(c4_angles.size())/angles_total << std::endl;
	  //---------------
	  c2_ratio.push_back( (double)(c2_angles.size())/angles_total );
	  c3_ratio.push_back( (double)(c3_angles.size())/angles_total );
	  c4_ratio.push_back( (double)(c4_angles.size())/angles_total );
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
  c3_angles_hist = makeRegularHistogram( c3_angles , hist_start , hist_end , nb_box );
  c4_angles_hist = makeRegularHistogram( c4_angles , hist_start , hist_end , nb_box );
  c2_ratio_hist = makeRegularHistogram( c2_ratio , min( c2_ratio ) , max( c2_ratio ) , 100 );
  c3_ratio_hist = makeRegularHistogram( c3_ratio , min( c3_ratio ) , max( c3_ratio ) , 100 );
  c4_ratio_hist = makeRegularHistogram( c4_ratio , min( c4_ratio ) , max( c4_ratio ) , 100 );
  //----------------------------------------------------
  // Normalizing for angles
  //----------------------------------------------------
  for ( int i=0 ; i < c2_angles_hist.size() ; i++ )
    {
      c2_angles_hist[i].value = (int)( c2_angles_hist[i].value/sin( center( c2_angles_hist[i] ) * M_PI/180  ) );
    }
  for ( int i=0 ; i < c3_angles_hist.size() ; i++ )
    {
      c3_angles_hist[i].value = (int)( c3_angles_hist[i].value/sin( center( c3_angles_hist[i] ) * M_PI/180 ) );
    }
  for ( int i=0 ; i < c4_angles_hist.size() ; i++ )
    {
      c4_angles_hist[i].value = (int)( c4_angles_hist[i].value/sin( center( c4_angles_hist[i] ) * M_PI/180 ) );
    }
  //----------------------------------------------------
  writeHistogram( c2_angles_out , normalizeHistogram( c2_angles_hist ) );
  writeHistogram( c3_angles_out , normalizeHistogram( c3_angles_hist ) );
  writeHistogram( c4_angles_out , normalizeHistogram( c4_angles_hist ) );
  writeHistogram( c2_ratio_out , normalizeHistogram( c2_ratio_hist ) );
  writeHistogram( c3_ratio_out , normalizeHistogram( c3_ratio_hist ) );
  writeHistogram( c4_ratio_out , normalizeHistogram( c4_ratio_hist ) );
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
  
  //--------------
  //Closing fluxes
  //----------------------
  input.close();
  c2_angles_out.close();
  c3_angles_out.close();
  c4_angles_out.close();
  c2_ratio_out.close();
  c3_ratio_out.close();
  c4_ratio_out.close();
  //----------------------
      
  return 0;
}
