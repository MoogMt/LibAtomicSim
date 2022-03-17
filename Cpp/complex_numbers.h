// Guard
#ifndef complex_number_h
#define complex_number_h

const double pi  = 3.1418;
const double tau = 2*pi;

// Complex Cartesian Form
//----------------------------------
class Complex
{
  protected:
    // Class variables
    //-------------------------
    double real;
    double imaginary;
    //-------------------------
  
  public:
    // Constructors
    Complex();
    Complex( double new_real );
    Complex( double new_real, double new_imaginary );
    // Accessors
    double re();
    double im();
    // Setters
    void re( double new_real );
    void im( double new_imaginary );
    void set( double new_real, double new_imaginary );
    // Conjugate
    Complex conj();
    // Module
    Complex mod( Complex complex );
};
// Operations on complex
Complex sum( Complex complex1, Complex complex2 );
Complex multiply( Complex complex1, Complex complex2 );
Complex* conj( Complex* complex_vector, int size );
//----------------------------------

// Complex Exponential Form
//----------------------------------
class ComplexExp
{
  protected:
    // Class variables
    //-------------------------
    double modulus;
    double phase;
    //-------------------------
  
  public:
    // Constructors
    ComplexExp();
    ComplexExp( double new_mod );
    ComplexExp( double new_mod, double new_phase );

    // Accessors
    double mod();
    double phi();

    // Setters
    void mod( double new_mod );
    void phi( double new_phase );
    void set( double new_mod, double new_phase );

    // Compute real and imaginary parts
    double re();
    double im();

    // Sums
    Complex sum( Complex complex1, Complex complex2 );

};
//----------------------------------



#endif
