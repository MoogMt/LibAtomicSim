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
  int start_step = 2000; // Start step
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
  
  //------------------------------------
  // Histograms
  //------------------------------------
  // Technical values
  double hist_start = 0.95;
  double hist_end   = 3.00;
  int nb_box = 400;
  //------------------------------------
  std::vector<Bin> hist_1CC; std::vector<double> CC1;
  std::vector<Bin> hist_2CC; std::vector<double> CC2;
  std::vector<Bin> hist_1CO; std::vector<double> CO1;
  std::vector<Bin> hist_2CO; std::vector<double> CO2;
  std::vector<Bin> hist_3CO; std::vector<double> CO3;
  std::vector<Bin> hist_4CO; std::vector<double> CO4;
  std::vector<Bin> hist_1OC; std::vector<double> OC1;
  std::vector<Bin> hist_2OC; std::vector<double> OC2;
  std::vector<Bin> hist_1OO; std::vector<double> OO1;
  //---------------------------------------------------

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
  std::ofstream outputCC_1nn ("1nearestCC.dat",  std::ios::out );
  std::ofstream outputCC_2nn ("2nearestCC.dat",  std::ios::out );
  std::ofstream outputCO_1nn ("1nearestCO.dat",  std::ios::out );
  std::ofstream outputCO_2nn ("2nearestCO.dat",  std::ios::out );
  std::ofstream outputCO_3nn ("3nearestCO.dat",  std::ios::out );
  std::ofstream outputCO_4nn ("4nearestCO.dat",  std::ios::out );  
  std::ofstream outputOC_1nn ("1nearestOC.dat",  std::ios::out );
  std::ofstream outputOC_2nn ("2nearestOC.dat",  std::ios::out );
  std::ofstream outputOO_1nn ("1nearestOO.dat",  std::ios::out );
  //--------------------------------------------------------------------------------

  //--------------
  // Atom Indexes
  //------------------------------
  std::vector<int> atom_indexesC;
  std::vector<int> atom_indexesO;
  //------------------------------

  std::cout <<"Starting..." << step << '\xd';

  //-------------------
  // Reading XYZ file
  //----------------------------------------------------
  while( readStepXYZfast( input , atom_list , lut_list, true, true ) )
    {
      if ( step % comp_step == 0 && step >= start_step && step <= end_step )
	{
	  if ( step == start_step )
	    {
	      atom_indexesC = getSpecieIndex( lut_list , "C" );
	      if ( atom_indexesC.size() == 0 )
		{
		  return 1;
		}
	      atom_indexesO = getSpecieIndex( lut_list , "O" );
	      if ( atom_indexesO.size() == 0 )
		{
		  return 1;
		}
	    }
	  makeContactMatrixDistance( cm , atom_list, cell , cut_off , lut_list );
	  appendVector( CC1 , getNNearest( cm , 1 , atom_indexesC , atom_indexesC ) );
	  appendVector( CC2 , getNNearest( cm , 2 , atom_indexesC , atom_indexesC ) );
	  appendVector( CO1 , getNNearest( cm , 1 , atom_indexesC , atom_indexesO ) );
	  appendVector( CO2 , getNNearest( cm , 2 , atom_indexesC , atom_indexesO ) );
	  appendVector( CO3 , getNNearest( cm , 3 , atom_indexesC , atom_indexesO ) );
	  appendVector( CO4 , getNNearest( cm , 4 , atom_indexesC , atom_indexesO ) );
	  appendVector( OC1 , getNNearest( cm , 1 , atom_indexesO , atom_indexesC ) );
	  appendVector( OC2 , getNNearest( cm , 2 , atom_indexesO , atom_indexesC ) );
	  appendVector( OO1 , getNNearest( cm , 1 , atom_indexesO , atom_indexesO ) );
	  std::cout << "Computing nearest neighbours, step: " << step << '\xd';
	}
      step++;
     }
  //----------------------------------------------------

  //----------------------
  // Writting histograms
  //--------------------------------------------------------------------------
  std::cout << "Making Histograms, this can take a while... " << std::endl;
  writeHistogram( outputCC_1nn , normalizeHistogram( makeRegularHistogram( CC1 , hist_start , hist_end , nb_box ) ) );
  writeHistogram( outputCC_2nn , normalizeHistogram( makeRegularHistogram( CC2 , hist_start , hist_end , nb_box ) ) );
  writeHistogram( outputCO_1nn , normalizeHistogram( makeRegularHistogram( CO1 , hist_start , hist_end , nb_box ) ) );
  writeHistogram( outputCO_2nn , normalizeHistogram( makeRegularHistogram( CO2 , hist_start , hist_end , nb_box ) ) );
  writeHistogram( outputCO_3nn , normalizeHistogram( makeRegularHistogram( CO3 , hist_start , hist_end , nb_box ) ) );
  writeHistogram( outputCO_4nn , normalizeHistogram( makeRegularHistogram( CO4 , hist_start , hist_end , nb_box ) ) );
  writeHistogram( outputOC_1nn , normalizeHistogram( makeRegularHistogram( OC1 , hist_start , hist_end , nb_box ) ) );
  writeHistogram( outputOC_2nn , normalizeHistogram( makeRegularHistogram( OC2 , hist_start , hist_end , nb_box ) ) );
  writeHistogram( outputOO_1nn , normalizeHistogram( makeRegularHistogram( OO1 , hist_start , hist_end , nb_box ) ) );
  //--------------------------------------------------------------------------

  //-----------------------------
  // Closing flux for histograms
  //-----------------------------------
  outputCC_1nn.close();
  outputCC_2nn.close();
  outputCO_1nn.close();
  outputCO_2nn.close();
  outputCO_3nn.close();
  outputCO_4nn.close();
  outputOC_1nn.close();
  outputOC_2nn.close();
  outputOO_1nn.close();
  //-----------------------------------
  
  //--------------
  //Closing fluxes
  //----------------------
  input.close();
   //----------------------

  std::cout << "All done! " << std::endl;
  
  return 0;
}
