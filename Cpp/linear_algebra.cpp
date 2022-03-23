#include <assert.h> 
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>

#include "complex_numbers.h"
#include "utils.h"
#include "linear_algebra.h"

// Basics
int computeIndex( int line, int col, int line_dim )
{
  return line*line_dim + col;
}

// Real
//------------------------------------------------------------------------------------------
// Init Matrix
Matrix::Matrix()
{
  size_x = 0;
  size_y = 0;
  elements = (double*) malloc( 0*sizeof(double) );
}
Matrix::Matrix( int size )
{
  size_x = size;
  size_y = size;
  elements=(double*) malloc( size*size*sizeof(double) );
  for( int i=0; i<size*size; i++ )
  {
    elements[i] = 0;
  }
}
Matrix::Matrix( int new_size_x, int new_size_y )
{
  size_x = new_size_x;
  size_y = new_size_y;
  elements=(double*) malloc( size_x*size_y*sizeof(double) );
  for( int i=0; i<size_x*size_y; i++ )
  {
    elements[i] = 0;
  }
}
// Accessors
int Matrix::sizex()
{
  return this->size_x;
}
int Matrix::sizey()
{
  return this->size_y;
}
double* Matrix::matrix()
{
  return this->elements;
}
double Matrix::element( int i, int j )
{
  return this->elements[ computeIndex(i,j,this->sizex())];
}
// Setters
void Matrix::sizex( int new_size_x )
{
  size_x = new_size_x;
  return;
}
void Matrix::sizey( int new_size_y )
{
  size_y = new_size_y;
  return;
}
void Matrix::matrix( double* matrix )
{
  int size_x1 = this->sizex();
  int size_y1 = this->sizey();
  for( int i=0; i<size_x1; i++ )
  {
    for( int j=0; j<size_y1; j++ )
    {
      this->elements[ computeIndex(i,j,size_x1 ) ] = matrix[ computeIndex(i,j,size_x1) ];
    }
  }
  return;
}
void Matrix::element( double element, int pos_x, int pos_y )
{
  int size_x1 = this->sizex();
  int size_y1 = this->sizey();
  this->elements[ computeIndex(pos_x,pos_y,size_x1 ) ] = element;
  return;
}
// Identity Matrix
void Matrix::eye( int size )
{
  size_x = size;
  size_y = size;
  elements=(double*) malloc( size*size*sizeof(double) );
  for( int i=0; i<size; i++ )
  {
    for( int j=0; j<size; j++ )
    {
      if( i == j )
      {
        elements[ computeIndex(i,j,size) ] = 1;
      }
      else
      {
        elements[ computeIndex(i,j,size) ] = 0;
      }
    }
  }
  return; 
}
std::ostream& operator<<( std::ostream& os, const Matrix& matrix )
{
  for( int j=0; j<matrix.size_y ; j++ )
  {
    for( int i=0; i<matrix.size_x ; i++ )
    {
      os << matrix.elements[ computeIndex(i,j,matrix.size_x)] << " ";
    }
    os << std::endl;
  }
  return os;
}
//==========================================================================================
double* zeroVector( int size )
{
  double* vector = (double*) malloc( size*sizeof(double) );
  for( int i=0; i<size; i++ )
  {
    vector[i] = 0;
  }
  return vector;
}
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
// Identity Matrix
double* eyeMatrix( int size )
{
  // Initialize matrix
  double* matrix = (double*) malloc( size*size*sizeof(double) );
  for( int i=0; i<size; i++ )
  {
    matrix[ computeIndex(i,i,size) ] = 1.0;
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
// Tranpose
double* transpose( double* matrix_origin, int size )
{
  double* matrix_tranpose = zeroMatrix( size );
  for( int i=0; i<size; i++ )
  {
    for( int j=0; j<size; j++ )
    {
      matrix_tranpose[ computeIndex(i,j,size) ] = matrix_origin[ computeIndex(j,i,size) ];
    }
  }
  return matrix_tranpose;
}
// Products
double  inner( double* vector1, double* vector2, int size )
{
  double inner=0;
  for( int i=0; i<size; i++ )
  {
    inner = inner +  vector1[i]*vector2[i];
  }
  return inner;
}
double* cross( double* vector1, double* vector2 )
{
  double* cross_vector = (double*) malloc( 3*sizeof(double) );
  cross_vector[0] = vector1[1]*vector2[2] - vector1[2]*vector2[1];
  cross_vector[1] = vector1[0]*vector2[2] - vector1[2]*vector2[0];
  cross_vector[2] = vector1[0]*vector2[1] - vector1[1]*vector2[0];
  return cross_vector;
} 
// Printing Matrix
void printMatrix( double* matrix, int size_x, int size_y )
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
//------------------------------------------------------------------------------------------

// Complex
//------------------------------------------------------------------------------------------
Complex* zeroVectorC( int size )
{
  Complex* vector = (Complex*) malloc( size*sizeof(Complex) );
  for( int i=0; i<size; i++ )
  {
    vector[i] = Complex(0,0);
  }
  return vector;
}
Complex* zeroMatrixC( int size_x, int size_y )
{
  int ndim=size_x*size_y;
  Complex* matrix = (Complex*) malloc( ndim*sizeof(Complex) );
  for( int i=0; i<ndim; i++ )
  {
    matrix[i] = Complex(0);
  }
  return matrix;
}
Complex* zeroMatrixC( int size )
{
  return zeroMatrixC( size, size);
}
Complex* eyeMatrixC( int size )
{
  Complex* matrix = (Complex*) malloc( size*size*sizeof(Complex) );
  for( int i=0; i<size; i++ )
  {
    matrix[ computeIndex(i,i,size) ] = Complex(1,0);
  }
  return matrix;
}
// Linear/Basic Operations
Complex* transposeC( Complex* matrix_origin, int size )
{
  Complex* matrix_transpose = zeroMatrixC( size );
  for( int i=0; i<size; i++ )
  {
    for( int j=0; j<size; j++ )
    {
      matrix_transpose[ computeIndex(i,j,size) ] = matrix_origin[ computeIndex(j,i,size) ];
    }
  }
  return matrix_transpose;
}
Complex* adjoint( Complex* matrix, int size )
{
  return transposeC( getConj( matrix, size ), size );
}
// Products
Complex* multiplyMatrix( Complex* matrix1, Complex* matrix2, int size1_x, int size_1y_2x, int size2_y)
{
  Complex* product_matrix = zeroMatrixC( size1_x, size2_y );
  for( int i=0; i<size1_x; i++ )
  {
    for( int j=0; j<size2_y; j++ )
    {
      for( int k=0; k<size_1y_2x; k++ )
      {
        product_matrix[ computeIndex(i,j,size1_x) ] = product_matrix[ computeIndex(i,j,size1_x) ] + matrix1[ computeIndex(i,k,size1_x) ]*matrix2[ computeIndex(k,j,size_1y_2x) ];
      }
    }
  }
  return product_matrix;
}
Complex inner( Complex* vector1, Complex* vector2, int size )
{
  Complex inner = Complex(0,0);
  for( int i=0; i<size; i++ )
  {
    inner = inner + vector1[i]*vector2[i];
  }
  return inner;
}
/*
Complex cross( Complex* vector1, Complex* vector2 )
{
  Complex* cross_vector = (Complex*) malloc( 3*sizeof(Complex) );
  cross_vector[0] = add( multiply( vector1[1], vector2[2] ), -multiply( vector1[2], vector2[1] ) );
  cross_vector[1] = add( multiply( vector1[0], vector2[2] ), -multiply( vector1[2], vector2[0] ) );
  cross_vector[2] = add( multiply( vector1[0], vector2[1] ), -multiply( vector1[1], vector2[0] ) );
  return cross_vector;
}*/
// Printing Matrix
void printMatrixC( Complex* matrix, int size_x, int size_y )
{
  for( int i=0; i<size_y; i++ )
  {
    int offset = size_x*i;
    for ( int j=0; j<size_x; j++ )
    {
      std::cout << matrix[ offset + j ].re();
      if( abs(matrix[ offset + j ].im()) > 0.00000000000001 )
      {
        std::cout << "+" << matrix[ offset + j ].im() << "i" << " ";      
      }
      else
      {
        std::cout << "   ";
      }
    }
    std::cout << std::endl;
  }
  return;
}
void printMatrixC( Complex* matrix, int size )
{
  printMatrixC( matrix, size, size);
  return;
}
//------------------------------------------------------------------------------------------