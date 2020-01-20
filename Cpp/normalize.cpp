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
  //----------
  // In File
  //-----------------------------------------
  std::string in_file_name;
  std::cout << "Name of the input file: " ;
  std::cin >> in_file_name;
  //-----------------------------------------
  
  //-----------------
  // Read histograms
  //---------------------------------------------------------------------------------------------
  std::vector<BinReal> hist = normalizeHistogram( readRegularHistogramReal( in_file_name ) );
  //---------------------------------------------------------------------------------------------

  //----------
  // Out File
  //---------------------------------------------------------------------------------------------
  std::string out_file_name;
  std::cout << "Name of the output file: " ;
  std::cin >> out_file_name;
  //---------------------------------------------------------------------------------------------
  
  //--------
  // Output
  //-------------------------------------
  std::ofstream output( out_file_name.c_str() );
  //-------------------------------------

  //----------------
  // Writting data
  //---------------------------------------------------------------------------------------------
  for ( int i=0 ; hist.size() ; i++ )
    {
      output << center( hist[i] ) << " " << hist[i].value << std::endl;
    }
  //---------------------------------------------------------------------------------------------

  //--------------
  //Closing fluxes
  //----------------------
  output.close();
  //----------------------
  
  return 0;
}
