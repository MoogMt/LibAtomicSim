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
    void set( Complex complex );
    // Conjugate  
    Complex getConj();
    void setConj();
    // Module
    double mod();
    
    // Self Operators
    Complex operator+( Complex complex2);
    Complex operator+( double real );
    void operator+=( Complex complex );
    void operator+=( double real );
    Complex operator-( Complex complex2);
    Complex operator-( double real );
    void operator-=( Complex complex );
    void operator-=( double real );
    Complex operator*( Complex complex2 );
    Complex operator*( double real );
    void operator*=( Complex complex );
    void operator*=( double real );
    Complex operator/( Complex complex2 );
    Complex operator/( double real );
    // Friend Operator
    friend std::ostream& operator << ( std::ostream& os, const Complex& complex );
    friend Complex operator-( Complex complex ); // Takng the opposite of a number
    friend Complex* getConj( Complex* complex_vector, int size ); 
    friend void setConj( Complex* complex_vector, int size );
};

// Complex Exponential Form
//----------------------------------
class ComplexExp
{
  protected:
    // Class variables
    //-------------------------
    double rho;
    double phi;
    //-------------------------
  
  public:
    // Constructors
    ComplexExp();
    ComplexExp( double new_mod );
    ComplexExp( double new_mod, double new_phase );

    // Accessors
    double mod();
    double phase();

    // Setters
    void mod( double new_mod );
    void phase( double new_phase );
    void set( double new_mod, double new_phase );

    // Compute real and imaginary parts
    double re();
    double im();
};
//----------------------------------



#endif
