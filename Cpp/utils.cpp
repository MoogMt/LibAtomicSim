#include <assert.h> 
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>

#include "complex_numbers.h"
#include "linear_algebra.h"
#include "utils.h"

using namespace std;

// Strings
//-------------------------------------------------------
std::string* initWordVector( int number )
{
  return (std::string*) malloc( number*sizeof( std::string ) );
}
//-------------------------------------------------------

// Error message
//-------------------------------------------------------
void indexOutOfBondError( int index, int max_index )
{
    std::cout << "Error!" << std::endl;
    std::cout << "Index out of bond: " << index << std::endl;
    std::cout << "Max index: " << max_index << std::endl;
}
//-------------------------------------------------------

// Cell Matrix
//-------------------------------------------------------
double* initCellMatrix()
{
  return eyeMatrix(3);
}
double* invertCellMatrix( double* matrix_to_invert )
{

  // Init inverse matrix
  double* inverse_matrix = (double*) malloc( 9*sizeof(double) );

  inverse_matrix[ computeIndex(0,0,3) ] = matrix_to_invert[ computeIndex(1,1,3) ]*matrix_to_invert[ computeIndex(2,2,3) ] \
  - matrix_to_invert[ computeIndex(2,1,3) ]*matrix_to_invert[ computeIndex(1,2,3) ];
  inverse_matrix[ computeIndex(1,0,3) ] = matrix_to_invert[ computeIndex(2,0,3) ]*matrix_to_invert[ computeIndex(1,0,3) ] \
  - matrix_to_invert[ computeIndex(1,0,3) ]*matrix_to_invert[ computeIndex(2,2,3) ];
  inverse_matrix[ computeIndex(2,0,3) ] = matrix_to_invert[ computeIndex(1,0,3) ]*matrix_to_invert[ computeIndex(2,1,3) ] \
  - matrix_to_invert[ computeIndex(2,0,3) ]*matrix_to_invert[ computeIndex(1,1,3) ];
  inverse_matrix[ computeIndex(0,1,3) ] = matrix_to_invert[ computeIndex(2,1,3) ]*matrix_to_invert[ computeIndex(0,2,3) ] \
  - matrix_to_invert[ computeIndex(0,1,3) ]*matrix_to_invert[ computeIndex(2,2,3) ];
  inverse_matrix[ computeIndex(1,1,3) ] = matrix_to_invert[ computeIndex(0,0,3) ]*matrix_to_invert[ computeIndex(2,2,3) ] \
  - matrix_to_invert[ computeIndex(2,0,3) ]*matrix_to_invert[ computeIndex(0,2,3) ];
  inverse_matrix[ computeIndex(2,1,3) ] = matrix_to_invert[ computeIndex(2,0,3) ]*matrix_to_invert[ computeIndex(0,1,3) ] \
  - matrix_to_invert[ computeIndex(0,0,3) ]*matrix_to_invert[ computeIndex(2,1,3) ];
  inverse_matrix[ computeIndex(0,2,3) ] = matrix_to_invert[ computeIndex(0,1,3) ]*matrix_to_invert[ computeIndex(1,2,3) ] \
  - matrix_to_invert[ computeIndex(1,1,3) ]*matrix_to_invert[ computeIndex(0,2,3) ];
  inverse_matrix[ computeIndex(1,2,3) ] = matrix_to_invert[ computeIndex(1,0,3) ]*matrix_to_invert[ computeIndex(0,2,3) ] \
  - matrix_to_invert[ computeIndex(0,0,3) ]*matrix_to_invert[ computeIndex(1,2,3) ];
  inverse_matrix[ computeIndex(2,2,3) ] = matrix_to_invert[ computeIndex(0,0,3) ]*matrix_to_invert[ computeIndex(1,1,3) ] \
  - matrix_to_invert[ computeIndex(1,2,3) ]*matrix_to_invert[ computeIndex(2,1,3) ];

  // Compute determinant
  double d = matrix_to_invert[ computeIndex(0,0,3) ]*inverse_matrix[ computeIndex(0,0,3) ]  \
  +   matrix_to_invert[ computeIndex(0,1,3) ]*inverse_matrix[ computeIndex(1,0,3) ]  \
  +   matrix_to_invert[ computeIndex(0,2,3) ]*inverse_matrix[ computeIndex(2,0,3) ];

  // Check that matrix is invertible (d!=0)
  try{ 
    // Precision is 10^-8
    if( abs(d) > 0.000000000001 )
    {
      multiplyby( inverse_matrix, 1.0/d, 3  );
      return inverse_matrix;
    }
    else
    {
      // Code of error is 0
      throw 0;
    }
  }
  // If not, throws error
  catch( int error )
  {
    std::cout << "Can't inverse matrix, determinant is null!";
    throw 0;
  }
}
//-------------------------------------------------------

// Geometry
//-------------------------------------------------------
double scalar( double* vector1, double* vector2, int ndim )
{
  double dist = 0;
  for ( int i=0; i<ndim; i++ )
  {
    dist = vector1[i]*vector2[i];
  }
  return dist;
}
//-------------------------------------------------------