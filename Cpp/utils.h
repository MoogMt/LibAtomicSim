// Guard
#ifndef utils_h
#define utils_h

using namespace std;

//-------------------------------------------------------
std::string* initWordVector( int number );
//-------------------------------------------------------


// Error message
//-------------------------------------------------------
void indexOutOfBondError( int index, int max_index );
//-------------------------------------------------------


// Matrix functions
//-------------------------------------------------------
// Basics
int computeIndex( int line, int col, int line_dim );
// Init matrixes
// -> Zero Matrix
double* zeroMatrix( int size_x, int size_y );
Complex* zeroMatrixC( int size_x, int size_y );
double* zeroMatrix( int size );
Complex* zeroMatrixC( int size );
// -> Identity Matrix
double* eyeMatrix( int size );
Complex* eyeMatrixC( int size );
// Linear Operations
void multiplyby( double* matrix, double scalar, int size_x, int size_y );
void multiplyby( double* matrix, double scalar, int size );
// Tranpose
double* transpose( double* matrix_origin, int size );
Complex* transpose( Complex* matrix_origin, int size );
// Printing
void printMatrix( double* matrix, int size_x, int size_y );
void printMatrix( double* matrix, int size_matrix );
void printMatrix( Complex* matrix, int size_x, int size_y );
//-------------------------------------------------------

// Cell Matrix
//-------------------------------------------------------
double* initCellMatrix();
double* invertCellMatrix( double* matrix_to_invert );
//-----------------------------------------------------------------------------

// Geometry
//-------------------------------------------------------
double scalar( double* vector1, double* vector2, int ndim );
//-------------------------------------------------------

#endif