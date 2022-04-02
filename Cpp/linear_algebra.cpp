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
  int ndim = size_x*size_y;
  double* matrix = (double*) malloc( ndim*sizeof(double) );
  for( int i=0; i<ndim; i++ )
  {
    matrix[i] = 0;
  }
  return matrix;
}
double* zeroMatrix( int size )
{
  return zeroMatrix( size, size );
}
// Identity Matrix
double* eyeMatrix( int size )
{
  // Initialize matrix
  double* matrix = (double*) malloc( size*size*sizeof(double) );
  for( int i=0; i<size; i++ )
  {
    matrix[ computeIndex( i, i, size ) ] = 1.0;
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
  return multiplyby( matrix, scalar, size, size );
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
        product_matrix[ computeIndex(i,j,size1_x) ] += matrix1[ computeIndex( i, k, size1_x ) ]*matrix2[ computeIndex( k, j, size_1y_2x ) ];
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
  return printMatrix( matrix, size_matrix, size_matrix );
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
Complex* eyeC( int size )
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
Complex* cross( Complex* vector1, Complex* vector2 )
{
  Complex* cross_vector = (Complex*) malloc( 3*sizeof(Complex) );
  cross_vector[0] = vector1[1]*vector2[2] - vector1[2]*vector2[1];
  cross_vector[1] = vector1[0]*vector2[2] - vector1[2]*vector2[0];
  cross_vector[2] = vector1[0]*vector2[1] - vector1[1]*vector2[0];
  return cross_vector;
}
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

// Real Vectors
//------------------------------------------------------------------------------------------
// Constructors
Vector::Vector()
{
  size = 0;
  direction = false;
  coefficients = (double*) malloc( 0*sizeof(double) );
}
Vector::Vector( bool new_dir )
{
  size = 0;
  direction = new_dir;
  coefficients = (double*) malloc( 0*sizeof(double) );
}
Vector::Vector( int new_size )
{
  size = new_size;
  coefficients = zeroVector( new_size );
  direction = false;
}
Vector::Vector( int new_size, bool new_dir )
{
  size = new_size;
  coefficients = zeroVector( new_size );
  direction = new_dir;
}
Vector::Vector( int new_size, double* new_coefficients )
{
  size = new_size;
  coefficients = (double*) malloc( new_size*sizeof(double) );
  for( int i=0 ; i<new_size ; i++ )
  {
    coefficients[i] = coefficients[i];
  }
  direction = false;
}
Vector::Vector( int new_size, double* new_coefficients, bool new_dir )
{
  size = new_size;
  coefficients = (double*) malloc( new_size*sizeof(double) );
  for( int i=0 ; i<new_size ; i++ )
  {
    coefficients[i] = new_coefficients[i];
  }
  direction = new_dir;
}
// Accessors
int Vector::getSize()
{
  return this->size;
}
double* Vector::vector()
{
  return this->coefficients;
}
double Vector::coef( int i )
{
  return this->coefficients[i];
}
bool Vector::dir()
{
  return this->direction;
}
// Setters
void Vector::setSize( int new_size )
{
  size = new_size;
  return;
}
void Vector::coef( double new_coef, int position )
{
  coefficients[ position ] = new_coef;
  return;
}
void Vector::vector( double* new_coefficients )
{
  for( int i=0; i<this->getSize(); i++ )
  {
    this->coef( new_coefficients[i], i );
  }
  return;
}
void Vector::dir( bool new_direction )
{
  direction=new_direction;
  return;
}
//------------------------------------------------------------------------------------------

// Complex Vectors
//------------------------------------------------------------------------------------------
// Constructors
VectorC::VectorC()
{
  size=0;
  direction = false;
  coefficients = (Complex*) malloc(0*sizeof(Complex) );
}
VectorC::VectorC( bool new_dir )
{
  size=0;
  direction = new_dir;
  coefficients = (Complex*) malloc(0*sizeof(Complex) );
}
VectorC::VectorC( int new_size )
{
  size = new_size;
  direction = false;
  coefficients = zeroVectorC( new_size );
}
VectorC::VectorC( int new_size, bool new_dir )
{
  size = new_size;
  direction = new_size;
  coefficients = zeroVectorC( new_size );
}
VectorC::VectorC( int new_size, Complex* new_coefficients)
{
  size = new_size;
  direction = false;  
  coefficients = new_coefficients;
}
VectorC::VectorC( int new_size, Complex* new_coefficients, bool new_dir )
{
  size      = new_size;
  direction = new_dir;  
  coefficients  = new_coefficients;
}
// Accessors
int VectorC::getSize()
{
  return this->size;
}
Complex* VectorC::vector()
{
  return this->coefficients;
}
Complex VectorC::coef( int i )
{
  return this->coefficients[i];
}
bool VectorC::dir()
{
  return this->direction;
}
// Setters
void VectorC::setSize( int new_size )
{
  size = new_size;
  return;
}
void VectorC::coef( Complex new_coef, int coef_position )
{
  coefficients[ coef_position ] = new_coef;
  return;
}
void VectorC::vector( Complex* new_coefficients )
{
  for( int i=0; i<this->getSize(); i++ )
  {
    this->coef( new_coefficients[i], i );
  }
  return;
}
void VectorC::dir( bool new_direction )
{
  direction = new_direction;
  return;
}
//------------------------------------------------------------------------------------------

// Real
//------------------------------------------------------------------------------------------
// Init Matrix
Matrix::Matrix()
{
  size_x = 0;
  size_y = 0;
  coefficients = (double*) malloc( 0*sizeof(double) );
}
Matrix::Matrix( int new_size )
{
  size_x = new_size;
  size_y = new_size;
  coefficients = zeroMatrix( new_size );
}
Matrix::Matrix( int new_size_x, int new_size_y )
{
  size_x = new_size_x;
  size_y = new_size_y;
  coefficients = zeroMatrix( new_size_x, new_size_y );
}
Matrix::Matrix( int new_size_x, int new_size_y, double* new_coefficients )
{
  size_x = new_size_x;
  size_y = new_size_y;
  coefficients=(double*) malloc( size_x*size_y*sizeof(double) );
  for( int i=0; i<size_x*size_y; i++ )
  {
    coefficients[i] = new_coefficients[i];
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
  return this->coefficients;
}
double Matrix::coef( int i, int j )
{
  return this->coefficients[ computeIndex( i, j, this->sizex() ) ];
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
void Matrix::coef( double new_coef, int coef_x, int coef_y )
{
  this->coef( new_coef, coef_x, coef_y );
  return;
}
void Matrix::matrix( double* matrix )
{
  for( int i=0; i<this->sizex(); i++ )
  {
    for( int j=0; j<this->sizey(); j++ )
    {
      this->coef( matrix[ computeIndex( i, j, this->sizex() ) ], i, j );
    }
  }
  return;
}
// Identity Matrix
void Matrix::eye( int size )
{
  size_x = size;
  size_y = size;
  coefficients = eyeMatrix( size );
  return; 
}
Matrix Matrix::operator+( Matrix matrix )
{
  Matrix new_matrix = Matrix( matrix.sizex() );
  try{
    if( sameDim( *this, matrix )  )
    {
      for( int i=0; i<matrix.sizex() ; i++ )
      {
        for( int j=0; j<matrix.sizey() ; j++ )
        {
          matrix.coef( this->coef(i,j) + matrix.coef(i,j), i, j );
        }
      }
    }
    else
    {
      throw(0);
    }
  }
  catch( int missmatch)
  {
    std::cout << "Wrong dimensions: " << std::endl;
    std::cout << "Matrix 1: " << this->sizex()   << " " << this->sizey() << std::endl;
    std::cout << "Matrix 2: " << matrix.sizex() << " " << matrix.sizey() << std::endl;
  }
  return new_matrix;
}
bool Matrix::isSymmetric()
{
  double eps=0.0000000000001;
  for( int i=1; i<this->sizex(); i++ )
  {
    for( int j=0; j<i-1; j++ )
    {
      if( abs( this->coef(i,j) - this->coef(j,i)) > eps )
      {
        return false;
      }
    }
  }
  return true;
}
// Friend Operators
bool sameDim( Matrix matrix1, Matrix matrix2 )
{
  if( matrix1.sizex() == matrix2.sizex() && matrix1.sizey() == matrix2.sizey() )
  {
    return true;
  }
  else return false;
}
std::ostream& operator<<( std::ostream& os, const Matrix& matrix )
{
  for( int j=0; j<matrix.size_y ; j++ )
  {
    for( int i=0; i<matrix.size_x ; i++ )
    {
      os << matrix.coefficients[ computeIndex(i,j,matrix.size_x)] << " ";
    }
    os << std::endl;
  }
  return os;
}
//------------------------------------------------------------------------------------------

// Complex Matrix
//------------------------------------------------------------------------------------------
// Constructors
MatrixC::MatrixC()
{
  size_x = 0;
  size_y = 0;
  coefficients = (Complex*) malloc( 0*sizeof(Complex) );
}
MatrixC::MatrixC( int new_size )
{
  size_x   = new_size;
  size_y   = new_size;
  coefficients = zeroVectorC( new_size );
}
MatrixC::MatrixC( int new_size, Complex* new_elements )
{
  size_x   = new_size;
  size_y   = new_size;
  coefficients = zeroMatrixC( new_size );
}
MatrixC::MatrixC( int new_size_x, int new_size_y )
{
  size_x = new_size_x;
  size_y = new_size_y;
  coefficients = zeroMatrixC( new_size_x, new_size_y );
}
MatrixC::MatrixC( int new_size_x, int new_size_y, Complex* new_elements )
{
  size_x = new_size_x;
  size_y = new_size_y;
  coefficients = zeroMatrixC( new_size_x, new_size_y );
}
// Accessors
int MatrixC::sizex()
{
  return this->size_x;
}
int MatrixC::sizey()
{
  return this->size_y;
}
Complex* MatrixC::matrix()
{
  return this->coefficients;
}
Complex MatrixC::coef( int i, int j )
{
  return this->coefficients[ computeIndex( i, j, this->sizex() ) ];
}
// Setters
void MatrixC::sizex( int new_size_x )
{
  size_x = new_size_x;
  return;
}
void MatrixC::sizey( int new_size_y )
{
  size_y = new_size_y;
  return;
}
void MatrixC::matrix( Complex* coefficients )
{
  int size_x1 = this->sizex();
  int size_y1 = this->sizey();
  for( int i=0; i<size_x1; i++ )
  {
    for( int j=0; j<size_y1; j++ )
    {
      this->coefficients[ computeIndex(i,j,size_x1 ) ] = coefficients[ computeIndex(i,j,size_x1) ];
    }
  }
  return;
}
void MatrixC::coef( Complex new_coef, int pos_x, int pos_y )
{
  int size_x1 = this->sizex();
  int size_y1 = this->sizey();
  this->coefficients[ computeIndex(pos_x,pos_y,size_x1 ) ] = new_coef;
  return;
}
// Identity Matrix
void MatrixC::eye( int size )
{
  size_x = size;
  size_y = size;
  coefficients = eyeC( size );
  return; 
}
MatrixC MatrixC::operator+( MatrixC matrix )
{
  MatrixC new_matrix = MatrixC( matrix.sizex() );
  try{
    if( sameDim( *this, matrix )  )
    {
      for( int i=0; i<matrix.sizex() ; i++ )
      {
        for( int j=0; j<matrix.sizey() ; j++ )
        {
          new_matrix.coef( this->coef(i,j) + matrix.coef(i,j), i, j );
        }
      }
    }
    else
    {
      throw( -1 );
    }
  }
  catch( int missmatch )
  {
    std::cout << "Wrong dimensions: " << std::endl;
    std::cout << "Matrix 1: " << this->sizex()   << " " << this->sizey() << std::endl;
    std::cout << "Matrix 2: " << matrix.sizex() << " " << matrix.sizey() << std::endl;
  }
  return new_matrix;
}
MatrixC MatrixC::operator-( MatrixC matrix )
{
  Complex* elements = (Complex*) malloc( sizeof( matrix.sizex()*matrix.sizey()*sizeof(Complex) ) );
  if( sameDim( *this, matrix )  )
  {
    for( int i=0; i<matrix.sizex() ; i++ )
    {
      for( int j=0; j<matrix.sizey() ; j++ )
      {
        coefficients[ computeIndex(i,j,matrix.sizex()) ] = this->coef(i,j) + matrix.coef(i,j);
      }
    }
    return MatrixC( matrix.sizex(), matrix.sizey(), elements );
  }
  else
  {
    return MatrixC( this->sizex(), this->sizey() );
  }
}
MatrixC MatrixC::transpose()
{
  MatrixC matrix_tranpose = MatrixC( this->sizex() );
  for( int i=0; i<this->sizex(); i++ )
  {
    for( int j=0; j<this->sizex(); j++ )
    {
      matrix_tranpose.coef( this->coef(j,i), i,j );
    }
  }
  return matrix_tranpose;
}
// Friend Operators
bool sameDim( MatrixC matrix1, MatrixC matrix2 )
{
  if( matrix1.sizex() == matrix2.sizex() && matrix1.sizey() == matrix2.sizey() )
  {
    return true;
  }
  else return false;
}
std::ostream& operator<<( std::ostream& os, const MatrixC& matrix )
{
  for( int j=0; j<matrix.size_y ; j++ )
  {
    for( int i=0; i<matrix.size_x ; i++ )
    {
      os << matrix.coefficients[ computeIndex(i,j,matrix.size_x)] << " ";
    }
    os << std::endl;
  }
  return os;
}
//------------------------------------------------------------------------------------------

