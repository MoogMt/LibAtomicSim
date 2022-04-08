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