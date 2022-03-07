#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>

#include "utils.h"
//#include "frame.h"
//#include "traj.h"

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
  
  double new_positions[3*3]{1, 2, 3, 4, 5, 6, 7, 8, 9};
  double new_positions2[2*3]{10,11,12,13,14,15};

  //matrix=(double*) malloc(9*sizeof(double) );
  //matrix_inv=(double*) malloc(9*sizeof(double) );

  //Frame frame(number);
  matrix=initCellMatrix();
  printMatrix(matrix,3);
  std::cout << "---------" << std::endl;
  matrix[ computeIndex(0,1,3) ] = 1;
  matrix[ computeIndex(0,2,3) ] = 1;
  matrix[ computeIndex(1,2,3) ] = 1;
  printMatrix(matrix,3);
  std::cout << "---------" << std::endl;
  matrix_inv=invertCellMatrix(matrix);
  printMatrix(matrix_inv,3);

  return 0;
}	       