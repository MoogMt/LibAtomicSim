#include <assert.h> 
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>

#include "complex_numbers.h"

#include "utils.h"

using namespace std;

// Strings
//-------------------------------------------------------
std::string* initWordVector( int number )
{
  return (std::string*) malloc( number*sizeof( std::string ) );
}
/*
std::string split( std::string str, std::string delim )
{

}*/
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

// Matrix functions
//-------------------------------------------------------
// Basics
int computeIndex( int line, int col, int line_dim )
{
  return line*line_dim + col;
}
// Init Matrix
double* zeroMatrix( int size_x, int size_y )
{
  int ndim=size_x*size_y;
  double* matrix = (double*) malloc( ndim*sizeof(double) );
  for( int i=0; i<ndim; i++ )
  {
    matrix[i] = 0;
  }

  return matrix;
}
double* zeroMatrix( int size )
{
  return zeroMatrix(size,size);
}
double* identityMatrix( int size )
{
  // Initialize matrix
  double* matrix = (double*) malloc( size*size*sizeof(double) );

  for( int i=0; i<size; i++ )
  {
    int offset=size*i;
    for( int j=0; j<size; j++ )
    {
      if(i==j)
      {
        matrix[offset+j]=1;
      }
      else
      {
        matrix[offset+j]=0;
      }
    }
  }
  return matrix;
}
// Linear Operations
void multiplyby( double* matrix, double scalar, int size_x, int size_y )
{
  int size_matrix = size_x*size_y;
  for( int element=0; element < size_matrix; element++ )
  {
    matrix[ element ] = scalar*matrix[ element ];
  }
  return;
}
void multiplyby( double* matrix, double scalar, int size )
{
  return multiplyby( matrix, scalar, size, size);
}
double* multiplyMatrix( double* matrix1, double* matrix2, int size1_x, int size_1y_2x, int size2_y)
{
  double* product_matrix = zeroMatrix( size1_x, size2_y );

  for( int i=0; i<size1_x; i++ )
  {
    for( int j=0; j<size2_y; j++ )
    {
      for( int k=0; k<size_1y_2x; k++ )
      {
        product_matrix[ computeIndex(i,j,size1_x) ] += matrix1[ computeIndex(i,k,size1_x) ]*matrix2[ computeIndex(k,j,size_1y_2x) ];
      }
    }
  }

  return product_matrix;
}
// Printing Matrix
void printMatrix( double* matrix, int size_y, int size_x )
{
  for( int i=0; i<size_y; i++ )
  {
    int offset = size_x*i;
    for ( int j=0; j<size_x; j++ )
    {
      std::cout << matrix[ offset + j ] << " ";      
    }
    std::cout << std::endl;
  }
  return;
}
void printMatrix( double* matrix, int size_matrix )
{
  return printMatrix(matrix,size_matrix,size_matrix);
}
//-------------------------------------------------------

// Cell Matrix
//-------------------------------------------------------
double* initCellMatrix()
{
  return identityMatrix(3);
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