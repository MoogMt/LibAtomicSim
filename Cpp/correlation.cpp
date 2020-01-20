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
  std::ifstream input( "TRAJEC.xyz" );
  //---------------------------------

  //--------
  // Output
  //--------------------------------------------------------------------------------
  std::ofstream corr( "corr.dat",  std::ios::out );
  std::ofstream check( "check.dat",  std::ios::out );
  //--------------------------------------------------------------------------------
  
  //----------------------
  // Physical parameters
  //---------------------------------------------------
  int step       = 1;      // Step counter
  int start_step = 5000;   // Start step
  int end_step   = 100000; // End Step 
  int comp_step  = 1;      // Frequency of computation
  int n_step = 0;
  int sim_stride = 5;     
  double timestep = 0.5;
  double femto=1e-15;
  double angstrom=1e-10;
  //---------------------------------------------------

  //------------------
  // Input from user
  //---------------------------------------------------
  std::cout << "Simulation Stride: ";
  std::cin >> sim_stride;
  double sim_timestep = sim_stride*timestep;
  std::cout << "Fraction to cut: ";
  double frac; 
  std::cin >> frac;
  //---------------------------------------------------
  
  //---------------
  // Initializers
  //--------------------------------------------------
  AtomList atom_list;      // Atoms in cell
  AtomList atom_list_prev; // Atoms in cell before
  std::vector<double> velocity_x; // Holding the velocity;
  std::vector<double> velocity_y; // Holding the velocity;
  std::vector<double> velocity_z; // Holding the velocity;
  AllTypeLUT lut_list; // LUT for types
  ContactMatrix cm;    // Contact Matrix
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

  //-------------------
  // Reading XYZ file
  //----------------------------------------------------
  while( readStepXYZfast( input , atom_list , lut_list, true, true ) )
    {
      if ( step % comp_step == 0 && step >= start_step && step <= end_step )
	{
	  if ( step != start_step )
	    {
	      for (int i=0 ; i < atom_list.x.size() ; i++ )
		{
		  velocity_x.push_back( ( atom_list.x[i] - atom_list_prev.x[i] )/sim_timestep*angstrom/femto );
		  velocity_y.push_back( ( atom_list.y[i] - atom_list_prev.y[i] )/sim_timestep*angstrom/femto );
		  velocity_z.push_back( ( atom_list.z[i] - atom_list_prev.z[i] )/sim_timestep*angstrom/femto );
		}
	    }
	  atom_list_prev = atom_list;
	  std::cout << "Computing velocities, step: " << step << std::endl;
	}
      step++;
     }
  //----------------------------------------------------

  //-----------------------
  // Computing Correlation
  //------------------------------------------------------------
   int nb_atoms = atom_list.x.size();
  n_step = velocity_x.size()/(double)(nb_atoms);
  std::vector<double> velocity_auto; velocity_auto.assign( n_step, 0. );
  std::vector<double> velocity_temp; velocity_temp.assign( n_step, 0. );
  // Compute autocorrelation velocity for each atom and averages the results
  for ( int i=0 ; i < nb_atoms ; i++)
    {
      std::cout << "Computing VACF for atom " << i << std::endl;
      // Calculation of the scalar product of velocities
      for ( int j=0 ; j < n_step ; j++ )
	{
	  double vx2 = velocity_x[ i ]*velocity_x[ i + j*nb_atoms ];
	  double vy2 = velocity_y[ i ]*velocity_y[ i + j*nb_atoms ];
	  double vz2 = velocity_z[ i ]*velocity_z[ i + j*nb_atoms ];
	  velocity_temp[j] =  vx2 + vy2 + vz2 ;
	}
      // Computation in itself
      autocorrelation( velocity_temp, frac );
      // Averages with other atoms
      for ( int j=0 ; j < n_step ; j++ )
	{
	  velocity_auto[j] += velocity_temp[j];
	}
    }
  // Normalize
  double max_auto = max(velocity_auto);
  for ( int i=0 ; i < velocity_auto.size() ; i++ )
    {
      velocity_auto[i] = velocity_auto[i]/max_auto;
    }
  //-------------------------------------------------------------- 

  //---------------
  // Writting data
  //--------------------------------------------------------------
  std::cout << "Writting VACF" << std::endl;
  for ( int i=0; i < velocity_auto.size() ; i++ )
    {
      corr << i*sim_timestep << " " << velocity_auto[i] << std::endl;
    }
  //---------------------------------------------------------------  
  
  //-----------------------------
  // Closing flux for histograms
  //-----------------------------------
  corr.close();
  check.close();
  //-----------------------------------
  
  //--------------
  //Closing fluxes
  //----------------------
  input.close();
   //----------------------
  
  return 0;
}
