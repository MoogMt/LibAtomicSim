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
double* zeroMatrix( int size_x, int size_y );
double* zeroMatrix( int size );
double* identityMatrix( int size );
// Linear Operations
void multiplyby( double* matrix, double scalar, int size_x, int size_y );
void multiplyby( double* matrix, double scalar, int size );
// Printing
void printMatrix( double* matrix, int size_y, int size_x );
void printMatrix( double* matrix, int size_matrix );
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