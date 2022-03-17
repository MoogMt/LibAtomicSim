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
// Printing
void printMatrix( double* matrix, int size_x, int size_y );
void printMatrix( double* matrix, int size_matrix );
//---------------------------------------------------------------

// Complex
//---------------------------------------------------------------
// Init 
Complex* zeroVector( int size );
Complex* zeroMatrixC( int size_x, int size_y );
Complex* zeroMatrixC( int size );
Complex* eyeMatrixC( int size );
// Basic Operatons
Complex* transpose( Complex* matrix_origin, int size );
Complex* adjoint( Complex* matrix, int size );
// Printing
void printMatrix( Complex* matrix, int size_x, int size_y );
void printMatrix( Complex* matrix, int size );
//---------------------------------------------------------------

#endif