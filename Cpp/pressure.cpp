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
#include "cell.h"
#include "histogram.h"
#include "sim.h"
//-------------------------

//================
// MAIN PROGRAM
//=====================================================================
int main( void )
{
  //--------
  // Input
  //---------------------------------
  std::ifstream press_in("STRESS");
  //--------------------------------

  //--------
  // Output
  //-------------------------------------
  std::ofstream press_out("pressure.dat");
  //-------------------------------------

  //---------------------------
  // Checking if STRESS exists
  //-------------------------------------------------------------
  if ( press_in )
    {
      std::cout << "File STRESS open" << std::endl;
    }
  else
    {
      std::cout << "File STRESS does not exists!" << std::endl;
      return 0;
    }
  //-------------------------------------------------------------
  
  //----------------------
  // Physical parameters
  //----------------------------------------------------
  int step       = 1;     // Step counter
  int start_step = 5000;  // Start step
  int end_step = 2000000; // End step
  int comp_step  = 1;     // Frequency of computation
  //----------------------------------------------------

  //--------------------------
  // Checking for user input
  //----------------------------------------------------------------------
  std::cout << "Give the starting step for data analysis: " std::endl;
  cin >> start_step;
  std::cout << "Give the end step for data analysis: " std::endl;
  cin >> end_step;
  std::cout << "Give the stride for data analysis: " std::endl;
  cin >> comp_step;
  //----------------------------------------------------------------------

  if ( start_step <= end_step )
    {
      std::cout << "OUT" << std::endl;
      return 1;
    }
  if ( comp_step > end_step-start_step )
    {
      std::cout << "OUT2" << std::endl;
      return 1;
    }

  //-------
  // Press 
  //----------------------------------------------------
  double pressure;
  std::vector<double> press_vec;
  std::vector<double> time_vec;
  //----------------------------------------------------

  //-------------------
  // Reading XYZ file
  //----------------------------------------------------
  if ( ! press_out )
    {
      //-------------
      // STRESS file
      //------------------------------------------------------------------------
      // Output file 
      std::ofstream press_time("pressure_time.dat");
      // Message to user
      std::cout << "Computing pressure" << step << std::endl;
      while(  readPressure( press_in, pressure ) )
	{
	  // Check over steps that are between start and end step
	  if ( step % comp_step == 0 && step >= start_step && step < end_step )
	    {
	      // Writting pressure as a function of time
	      press_time << step - start_step << " " <<  pressure << std::endl;
	      // Stocking pressure value in memory
	      press_vec.push_back( pressure );
	      // Message for user
	      std::cout << "Step = " << step << '\xd';
	    }      
	  step++;
	}
      std::cout << "Done" << std::endl;
      //------------------------------------------------------------------------

      //---------------------------------------
      // Making histogram and writting to file
      //------------------------------------------------------------------------
      writeHistogram( press_out , normalizeHistogram( makeRegularHistogram( press_vec , min( press_vec ) , max( press_vec) , 200 ) ) );
      //------------------------------------------------------------------------
      
      press_time.close();
    }
  else
    {
      // Output file 
      std::ifstream press_time("pressure_time.dat");
      //---------------------------
      // Reading pressure.dat file
      //-----------------------------------------------------
      std::cout << "Reading pressure_time.dat" << std::endl;
      while ( readPressure( press_time , pressure ) )
	{
	  press_vec.push_back( pressure);
	}
      //-----------------------------------------------------
      std::cout << "Done!" << std::endl;
      press_time.close();
    }
  //---------------------------------------------------------

  //----------------------------------------------
  // Print average and variance of the pressure 
  //---------------------------------------------------------------
  std::cout << "Average Pressure: " << average( press_vec )  << std::endl;
  //---------------------------------------------------------------
  
  //--------------
  //Closing fluxes
  //----------------------
  press_in.close();
  press_out.close();
  //----------------------
  
  return 0;
}
