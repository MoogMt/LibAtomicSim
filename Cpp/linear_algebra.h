// Guard
#ifndef lin_al
#define lin_al

using namespace std;

// Basics
//---------------------------------------------------------------
int computeIndex( int line, int col, int line_dim );
//---------------------------------------------------------------

//----------------------------------
class Vector
{
  protected: 
    // Variables
    double* elements;
    int size;
    bool direction;
  public:
    // Constructors
    Vector();
    Vector( bool new_dir );
    Vector( int new_size );
    Vector( int new_size, bool new_dir );
    Vector( int new_size, double* new_elements );
    Vector( int new_size, double* new_elements, bool new_dir );
    // Accessors
    int getSize();
    double* vector();
    double element( int i );
    bool dir();
    // Setters
    void setSize( int size );
    void element( double element, int position );
    void vector( double* element );
    void dir( bool new_direction );
};
class VectorC
{
  protected:
    // Variables
    Complex* elements;
    int  size;
    bool dir;
};
class Matrix
{
  protected:
    // Class variables
    double* elements;
    int size_x;
    int size_y;
  public:
    // Constructors
    Matrix();
    Matrix( int size );
    Matrix( int size_x, int size_y );
    Matrix( int new_size_x, int new_size_y, double* new_elements );
    // Accessors
    int sizex();
    int sizey();
    double* matrix();
    double element( int i, int j );
    // Setters
    void sizex( int size_x );
    void sizey( int size_y );
    void matrix( double* matrix );
    void element( double element, int pos_x, int pos_y );

    // Self Operators
    // Init
    void eye( int size );
    // Basic Operators
    Matrix operator+( Matrix matrix );
    // Friend Operator
    friend std::ostream& operator<<( std::ostream& os, const Matrix& matrix );
    friend bool sameDim( Matrix matrix1, Matrix matrix2 );
};
class MatrixC
{
  public:
    // Class variables
    Complex* elements;
    int size_x;
    int size_y;
};
//----------------------------------

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
double* cross( double* vector1, double* vector2 );
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