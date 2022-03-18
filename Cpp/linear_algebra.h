// Guard
#ifndef lin_al
#define lin_al

using namespace std;

// Basics
//---------------------------------------------------------------
int computeIndex( int line, int col, int line_dim );
//---------------------------------------------------------------

// Real
//---------------------------------------------------------------
// -> Zero Matrix
double* zeroVector( int size );
double* zeroMatrix( int size_x, int size_y );
double* zeroMatrix( int size );
// -> Identity Matrix
double* eyeMatrix( int size );
// Linear Operations
void multiplyby( double* matrix, double scalar, int size_x, int size_y );
void multiplyby( double* matrix, double scalar, int size );
// Tranpose
double* transpose( double* matrix_origin, int size );
// Basic Operations
double inner( double* vector1, double* vector2, int size );
double cross( double* vector1, double* vector2 );
// Printing
void printMatrix( double* matrix, int size_x, int size_y );
void printMatrix( double* matrix, int size_matrix );
//---------------------------------------------------------------

// Complex
//---------------------------------------------------------------
// Init 
Complex* zeroVectorC( int size );
Complex* zeroMatrixC( int size_x, int size_y );
Complex* zeroMatrixC( int size );
Complex* eyeMatrixC( int size );
// Basic Operatons
Complex* transpose( Complex* matrix_origin, int size );
Complex* adjoint( Complex* matrix, int size );
// Products
Complex* multiplyMatrix( Complex* matrix1, Complex* matrix2, int size1_x, int size_1y_2x, int size2_y);
Complex inner( Complex* vector1, Complex* vector2, int size );
// Printing
void printMatrix( Complex* matrix, int size_x, int size_y );
void printMatrix( Complex* matrix, int size );
//---------------------------------------------------------------

#endif