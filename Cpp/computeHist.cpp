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
  std::ofstream output("pv.dat");
  std::ofstream output2("var.dat");
  //-------------------------------------
  
  //-----------------
  // Read histograms
  //---------------------------------------------------------------------------------------------
  std::vector<BinReal> press_882_2000_hist = readRegularHistogramReal("8.82/2000K/pressure.dat");
  std::vector<BinReal> press_900_2000_hist = readRegularHistogramReal("9.0/2000K/pressure.dat");
  std::vector<BinReal> press_910_2000_hist = readRegularHistogramReal("9.1/2000K/pressure.dat");
  std::vector<BinReal> press_980_2000_hist = readRegularHistogramReal("9.8/2000K/pressure.dat");
  std::vector<BinReal> press_882_2500_hist = readRegularHistogramReal("8.82/2500K/pressure.dat");
  std::vector<BinReal> press_900_2500_hist = readRegularHistogramReal("9.0/2500K/pressure.dat");
  std::vector<BinReal> press_980_2500_hist = readRegularHistogramReal("9.8/2500K/pressure.dat");
  std::vector<BinReal> press_882_3000_hist = readRegularHistogramReal("8.82/3000K/pressure.dat");
  std::vector<BinReal> press_900_3000_hist = readRegularHistogramReal("9.0/3000K/pressure.dat");
  std::vector<BinReal> press_980_3000_hist = readRegularHistogramReal("9.8/3000K/pressure.dat");
  //---------------------------------------------------------------------------------------------

  //------------------
  // Compute averages
  //---------------------
  output << 8.82*8.82*8.82/96 << " " << average(press_882_2000_hist) << " " << average(press_882_2500_hist) << " " << average(press_882_3000_hist) << " " << variance(press_882_2000_hist) << " " << variance(press_882_2500_hist) << " " << variance(press_882_3000_hist) << std::endl;
  output << 9.0*9*9/96 << " " << average(press_900_2000_hist) << " " << average(press_900_2500_hist) << " " << average(press_900_3000_hist) << " " << variance(press_900_2000_hist) << " " << variance(press_900_2500_hist) << " " << variance(press_900_3000_hist) <<std::endl;
  output << 9.8*9.8*9.8/96 << " " << average(press_980_2000_hist) << " " << average(press_980_2500_hist) << " " << average(press_980_3000_hist) << " " <<  variance(press_980_2000_hist) << " " << variance(press_980_2500_hist) << " " << variance(press_980_3000_hist) << std::endl;
  //---------------------
  
  //--------------
  //Closing fluxes
  //----------------------
  input.close();
  output.close();
  output2.close();
  //----------------------
  
  return 0;
}
