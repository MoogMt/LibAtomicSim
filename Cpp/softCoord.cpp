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

  //----------------------
  // Physical parameters
  //--------------------------------------
  int step       = 1;  // Step counter
  int start_step = 5000; // Start step
  int end_step   = 10000000; // End Step 
  int comp_step  = 1; // Frequency of computation
  //--------------------------------------

  //---------------
  // Initializers
  //--------------------------------------------------
  AtomList  atom_list;  // Atoms in cell
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

  //-----------------
  // Reading Cut-Off
  //-------------------------------------------------------------------
  CutOffMatrix cut_off;
  if ( ! readCutOff( "cut_off.dat" , cut_off , lut_list ) )
    {
      return 1;
    }
  //-------------------------------------------------------------------
  
  //--------
  // Output
  //--------------------------------------------------------------------------------
  std::ofstream matrixC_out("C_coord_hist.dat",  std::ios::out );
  std::ofstream matrixO_out("O_coord_hist.dat",  std::ios::out );
  std::ofstream matrixAll_out("All_coord_hist.dat",  std::ios::out );
  std::ofstream matrixAll_stock("All_coord.dat",  std::ios::out );
  //--------------------------------------------------------------------------------

  //--------------
  // Atom Indexes
  //-------------------------
  std::vector<int> atom_indexesC;
  std::vector<int> atom_indexesO;
  //-------------------------

  std::vector<double> coord_C;
  std::vector<double> coord_O;
  std::vector<double> coord_total;
  double sum=0;
  
  //-------------------
  // Reading XYZ file
  //----------------------------------------------------
  while( readStepXYZfast( input , atom_list , lut_list, true, true ) )
    {
      if ( step % comp_step == 0 && step >= start_step && step <= end_step )
	{
	  cm  = makeContactMatrixSoft( atom_list , cell , cut_off , lut_list , 1.75 , 10 , 50 );
	  if ( step == start_step )
	    {
	      atom_indexesC = getSpecieIndex( lut_list, "C" );
	      atom_indexesO = getSpecieIndex( lut_list, "O" );
	    }
	  matrixAll_stock << step << " ";
	  for ( int i = 0 ; i < atom_indexesC.size() ; i++ )
	    {
	      sum = cumSum( getAtomContact( cm , atom_indexesC[i] ) );
	      coord_C.push_back( sum );
	      coord_total.push_back( sum );
	      matrixAll_stock << sum << " ";
	    }
	  for ( int i = 0 ; i < atom_indexesO.size() ; i++ )
	    {
	      sum = cumSum( getAtomContact( cm , atom_indexesO[i] ) );
	      coord_O.push_back( sum );
	      coord_total.push_back( sum );
	      matrixAll_stock << sum << " ";
	    }
	  matrixAll_stock << std::endl;
	  std::cout << "Reading file, step: "<< step << std::endl;
	}
      step++;
    }
  //----------------------------------------------------

  //-------------------
  // MAKING HISTOGRAMS
  //----------------------------------------------------
  int nb_box = 200;
  std::cout << "Making histograms: This can take a while...." << std::endl;
  writeHistogram( matrixC_out , normalizeHistogram( makeRegularHistogram( coord_C , min( coord_C ) , max( coord_C ) , nb_box ) ) );
  writeHistogram( matrixO_out , normalizeHistogram( makeRegularHistogram( coord_O , min( coord_O ) , max( coord_O ) , nb_box ) ) );
  writeHistogram( matrixAll_out , normalizeHistogram( makeRegularHistogram( coord_total , min( coord_total ) , max( coord_total ) , nb_box ) ) );
  //----------------------------------------------------
  
  //---------------------
  // Closing output flux
  //---------------------
  matrixC_out.close();
  matrixO_out.close();
  matrixAll_out.close();
  //-----------------------------------
  
  //--------------
  //Closing fluxes
  //----------------------
  input.close();
  //----------------------

  return 0;
}
