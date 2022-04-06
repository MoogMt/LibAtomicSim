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
double* zeroMatrix( int new_size );
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

//----------------------------------
class Vector
{
  protected: 
    // Variables
    int size;
    bool direction;
    double* coefficients;
  public:
    // Constructors
    Vector();
    Vector( bool new_dir );
    Vector( int new_size );
    Vector( int new_size, bool new_dir );
    Vector( int new_size, double* new_coefs );
    Vector( int new_size, double* new_coefs, bool new_dir );
    // Accessors
    int getSize();
    double* vector();
    double coef( int i );
    bool dir();
    // Setters
    void setSize( int new_size );
    void coef( double new_coef, int position );
    void vector( double* new_vector );
    void dir( bool new_direction );
    // Operator
    Vector operator+( Vector vector1 );
    Vector operator-( Vector vector1 );
    Vector operator*( Vector vector1 );
};
class VectorC
{
  protected:
    // Variables
    int  size;
    bool direction;
    Complex* coefficients;
  public:
    // Constructors
    VectorC();
    VectorC( bool new_dir );
    VectorC( int new_size );
    VectorC( int new_size, bool new_dir );
    VectorC( int new_size, Complex* new_coefs );
    VectorC( int new_size, Complex* new_coefs, bool new_dir );
    // Accessors
    int getSize();
    Complex* vector();
    Complex coef( int i );
    bool dir();
    // Setters
    void setSize( int size );
    void coef( Complex new_coef, int coef_pos );
    void vector( Complex* new_coef );
    void dir( bool new_direction );
};
class Matrix
{
  protected:
    // Class variables
    int size_x;
    int size_y;
    double* coefficients;
  public:
    // Constructors
    Matrix();
    Matrix( int size );
    Matrix( int size_x, int size_y );
    Matrix( int new_size_x, int new_size_y, double* new_coefficients );
    // Accessors
    int sizex();
    int sizey();
    double* matrix();
    double* row( int row_nb );
    double* col( int col_nb );
    double coef( int i, int j );
    // Setters
    void sizex( int size_x );
    void sizey( int size_y );
    void matrix( double* matrix );
    void row( int row_nb, double* new_row );
    void col( int col_nb, double* new_col );
    void coef( double new_coef, int coef_x, int coef_y );
    // Self Operators
    // Init
    void eye( int size );
    // Basic Operators
    Matrix operator+( Matrix matrix );
    bool isSymmetric();
    // Friend Operator
    friend std::ostream& operator<<( std::ostream& os, const Matrix& matrix );
    friend bool sameDim( Matrix matrix1, Matrix matrix2 );
};
class MatrixC
{
  protected:
    // Class variables
    int size_x;
    int size_y;
    Complex* coefficients;
  public:
    // Constructors
    MatrixC();
    MatrixC( int size );
    MatrixC( int size, Complex* new_coefficients );
    MatrixC( int size_x, int size_y );
    MatrixC( int new_size_x, int new_size_y, Complex* new_coefficients );
    // Accessors    void row( int row_nb, double* new_row );
    void col( int col_nb, double* new_col );
    int sizex();
    int sizey();
    Complex* matrix();
    Complex* row( int row_nb );
    Complex* col( int col_nb );
    Complex coef( int i, int j );
    // Setters
    void sizex( int size_x );
    void sizey( int size_y );
    void matrix( Complex* matrix );
    void row( int row_nb, Complex* new_row );
    void col( int col_nb, Complex* new_col );
    void coef( Complex new_coef, int coef_x, int coef_y );
    // Self Operators
    // Init
    void eye( int size );
    // Basic Operators
    MatrixC operator+( MatrixC matrix );
    MatrixC operator-( MatrixC matrix );
    MatrixC transpose();
    // Friend Operator
    friend std::ostream& operator<<( std::ostream& os, const MatrixC& matrix );
    friend bool sameDim( MatrixC matrix1, MatrixC matrix2 );
};
//----------------------------------

#endif