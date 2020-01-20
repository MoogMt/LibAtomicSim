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
  std::ifstream input("mo4s8.xyz");
  //--------------------------------

  //--------
  // Output
  //-------------------------------------------------
  std::ofstream output_sprint("distance_sprint.dat");
  std::ofstream output_real("distance_real.dat");
  std::ofstream output_sprint_hist("distance_histo.dat");
  //-------------------------------------------------
  
  //----------------------
  // Physical parameters
  //--------------------------------------
  int step       = 1;  // Step counter
  int start_step = 1; // Start step
  int end_step   = 140000;
  int comp_step  = 1; // Frequency of computation
  //--------------------------------------

  //---------------
  // Initializers
  //--------------------------------------------------------
  AtomList  atom_list_ini;          // Atoms in cell
  AtomList  atom_list;              // Atoms in cell
  AllTypeLUT lut_list;              // LUT for types
  ContactMatrix cm_ini;             // Contact Matrix Ini
  ContactMatrix cm_i;               // Contact Matrix at i
  std::vector<double> distance_vec; // Distance vector
  //--------------------------------------------------------

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

  //------------------
  // SPRINT variables
  //----------------------------------
  double r0 = 2.3812976205;
  int n = 8, m = 24;
  //----------------------------------
  
  //-------------------
  // Reading XYZ file
  //----------------------------------------------------
  std::cout << "Computing SPRINT distance: " << std::endl;
  // Loop over XYZ steps
  while( readStepXYZfast( input , atom_list , lut_list, true, true ) )
    {
      if ( step == start_step )
	{
	  atom_list_ini = atom_list;
	  cm_ini = makeContactMatrixSoft ( atom_list , cell , cut_off , lut_list , r0 , n,  m );
	}
      else if ( step % comp_step == 0 && step > start_step && step <= end_step )
	{
	  // Makes the contact matrix for step i
	  cm_i = makeContactMatrixSoft ( atom_list , cell , cut_off , lut_list , r0 , n,  m );
	  // Computes difference between contact matrix
	  cm_i.matrix = difference( cm_i.matrix , cm_ini.matrix );
	  double distance = norm(cm_i.matrix);
	  // Prints norm of resulting difference
	  output_sprint << step << " " << distance<< std::endl;
	  // Adding to vector
	  distance_vec.push_back(distance);
	  // Compute distance real
	  distance=0;
	  for ( int i=0 ; i < atom_list.x.size() ; i++ )
	    {
	      double x = atom_list.x[i] - atom_list_ini.x[i];
	      double y = atom_list.y[i] - atom_list_ini.y[i];
	      double z = atom_list.z[i] - atom_list_ini.z[i];
	      distance += sqrt( x*x + y*y + z*z );
	    }
	  // Normalize
	  distance = distance/(double)(atom_list.x.size());
	  // Writting distance
	  output_real << step << " " << distance << std::endl;
	  // Message out
	  std::cout << step << "\xd";
	}
      else if ( step >= end_step ) break;
      //------------------------------------------------
      step++;
     }
  std::cout << "Done !" << std::endl;
  //----------------------------------------------------

  //--------------------
  // Making Histogram
  //---------------------------------------------------------------------------------
  writeHistogram( output_sprint_hist , normalizeHistogram( makeRegularHistogram( distance_vec , min(distance_vec) , max(distance_vec) , 200 ) ) );
  //---------------------------------------------------------------------------------
  
  //--------------
  //Closing fluxes
  //----------------------
  input.close();
  output_sprint.close();
  output_real.close();
  output_sprint_hist.close();
  //----------------------
      
  return 0;
}
