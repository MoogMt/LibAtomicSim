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
  std::ifstream input("histogram.dat");
  //--------------------------------

  //--------
  // Output
  //-------------------------------------
  std::ofstream output2("angles2.dat");
  std::ofstream output3("angles3.dat");
  std::ofstream output4("angles4.dat");
  //-------------------------------------
  
  //-----------------
  // Read histograms
  //---------------------------------------------------------------------------------------------
  // CO2
  //---------------------------------------------------------------------------------------------
  std::vector<BinReal> angles2_882_2000_hist = readRegularHistogramReal("8.82/2000K/c2angles.dat");
  std::vector<BinReal> angles2_900_2000_hist = readRegularHistogramReal("9.0/2000K/c2angles.dat");
  std::vector<BinReal> angles2_980_2000_hist = readRegularHistogramReal("9.8/2000K/c2angles.dat");
  std::vector<BinReal> angles2_882_2500_hist = readRegularHistogramReal("8.82/2500K/c2angles.dat");
  std::vector<BinReal> angles2_900_2500_hist = readRegularHistogramReal("9.0/2500K/c2angles.dat");
  std::vector<BinReal> angles2_980_2500_hist = readRegularHistogramReal("9.8/2500K/c2angles.dat");
  std::vector<BinReal> angles2_882_3000_hist = readRegularHistogramReal("8.82/3000K/c2angles.dat");
  std::vector<BinReal> angles2_900_3000_hist = readRegularHistogramReal("9.0/3000K/c2angles.dat");
  std::vector<BinReal> angles2_980_3000_hist = readRegularHistogramReal("9.8/3000K/c2angles.dat");
  //---------------------------------------------------------------------------------------------
  // CO3
  //---------------------------------------------------------------------------------------------
  std::vector<BinReal> angles3_882_2000_hist = readRegularHistogramReal("8.82/2000K/c3angles.dat");
  std::vector<BinReal> angles3_900_2000_hist = readRegularHistogramReal("9.0/2000K/c3angles.dat");
  std::vector<BinReal> angles3_882_2500_hist = readRegularHistogramReal("8.82/2500K/c3angles.dat");
  std::vector<BinReal> angles3_900_2500_hist = readRegularHistogramReal("9.0/2500K/c3angles.dat");
  std::vector<BinReal> angles3_882_3000_hist = readRegularHistogramReal("8.82/3000K/c3angles.dat");
  std::vector<BinReal> angles3_900_3000_hist = readRegularHistogramReal("9.0/3000K/c3angles.dat");
  //---------------------------------------------------------------------------------------------
  // CO4
  //---------------------------------------------------------------------------------------------
  std::vector<BinReal> angles4_882_2000_hist = readRegularHistogramReal("8.82/2000K/c4angles.dat");
  std::vector<BinReal> angles4_900_2000_hist = readRegularHistogramReal("9.0/2000K/c4angles.dat");
  std::vector<BinReal> angles4_882_2500_hist = readRegularHistogramReal("8.82/2500K/c4angles.dat");
  std::vector<BinReal> angles4_900_2500_hist = readRegularHistogramReal("9.0/2500K/c4angles.dat");
  std::vector<BinReal> angles4_882_3000_hist = readRegularHistogramReal("8.82/3000K/c4angles.dat");
  std::vector<BinReal> angles4_900_3000_hist = readRegularHistogramReal("9.0/3000K/c4angles.dat");
  //---------------------------------------------------------------------------------------------

  //------------------
  // Compute averages
  //---------------------
  output2 << "8.82 " << average( angles2_882_2000_hist ) << " " << average( angles2_882_2500_hist ) << " " << average( angles2_882_3000_hist ) << " " << std::endl;
  output2 << "9.00 " << average( angles2_900_2000_hist ) << " " << average( angles2_900_2500_hist ) << " " << average( angles2_900_3000_hist ) << " " << std::endl;
  output2 << "9.80 " << average( angles2_980_2000_hist ) << " " << average( angles2_980_2500_hist ) << " " << average( angles2_980_3000_hist ) << " " << std::endl;
  //---------------------
    output3 << "8.82 " << average( angles3_882_2000_hist ) << " " << average( angles3_882_2500_hist ) << " " << average( angles3_882_3000_hist ) << " " << std::endl;
    output3 << "9.00 " << average( angles3_900_2000_hist ) << " " << average( angles3_900_2500_hist ) << " " << average( angles3_900_3000_hist ) << " " << std::endl;
  //------------------------
    output4 << "8.82 " << average( angles4_882_2000_hist ) << " " << average( angles4_882_2500_hist ) << " " << average( angles4_882_3000_hist ) << " " << std::endl;
    output4 << "9.00 " << average( angles4_900_2000_hist ) << " " << average( angles4_900_2500_hist ) << " " << average( angles4_900_3000_hist ) << " " << std::endl;
  //------------------------
    
  //--------------
  //Closing fluxes
  //----------------------
  input.close();
  output2.close();
  output3.close();
  output4.close();
  //----------------------
  
  return 0;
}
