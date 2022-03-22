#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>

#include "general_io.h"
#include "complex_numbers.h"
#include "linear_algebra.h"
#include "utils.h"
#include "frame.h"
#include "traj.h"

using namespace std;

// Main Program
//-------------------------------------------------------
int main()
{
  int number=3;
  int number2=2;  

  double* matrix;
  double* matrix_inv;

  std::string new_names[3]{"C","O","O"};
  std::string new_names2[2]{"C2","O2"} ;
  std::string file_path;

  bool status=false;
  
  file_path="/home/moogmt/file_test.dat";

  int number_line=0;
  double new_positions2[2*3]{10,11,12,13,14,15};

  //matrix=(double*) malloc(9*sizeof(double) );
  //matrix_inv=(double*) malloc(9*sizeof(double) );

  //Frame frame(number);
  matrix=initCellMatrix();
  printMatrix(matrix,3);
  matrix[ computeIndex(0,1,3) ] = 1;
  matrix[ computeIndex(0,2,3) ] = 1;
  matrix[ computeIndex(1,2,3) ] = 1;
  printMatrix(matrix,3);
  std::cout << "---------" << std::endl;
  matrix_inv=invertCellMatrix(matrix);
  printMatrix(matrix_inv,3);

  status = fileExists( file_path );
  if( status )
  {
    std::cout << "yay" << std::endl;
    number_line = getNumberLine( file_path );
    std::cout << "File has " << number_line << " lines." << std::endl;
  }
  else 
  {
    std::cout << "nay :(" << std::endl;
  }

  Complex test(0,1);

  std::cout << test << std::endl;

  std::cout << M_PI << std::endl;

  std::cout << std::fmod(6.0,M_PI) << std::endl;
  return 0;
}	       