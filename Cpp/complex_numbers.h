// Guard
#ifndef complex_number_h
#define complex_number_h

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
    double getReal();
    double getIm();
    // Setters
    void setReal( double new_real );
    void setIm( double new_imaginary );
    void setNumber( double new_real, double new_imaginary );

    // Compute conjugate
    Complex getConjugate();
};
// Operations on complex
Complex sum( Complex complex1, Complex complex2 );
Complex multiply( Complex complex1, Complex complex2 );
Complex* getConjugate( Complex* complex_vector, int size );
//----------------------------------

// Complex Exponential Form
//----------------------------------
class ComplexExp
{
  protected:
    // Class variables
    //-------------------------
    double mods;
    double phase;
    //-------------------------
  
  public:
    // Constructors
    ComplexExp();
    ComplexExp( double new_modulus );
    ComplexExp( double new_modulus, double new_phase );

    // Accessors
    double getModulus();
    double getPhase();

    // Setters
    void setModulus( double new_mod );
    void setPhase( double new_phase );
    void setNumber( double new_mod, double new_phase );

    Complex sum( Complex complex1, Complex complex2 );
};
//----------------------------------

#endif
