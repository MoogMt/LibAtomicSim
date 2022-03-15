#include <assert.h> 
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>

#include "complex_numbers.h"

// Complex cartesian form
//----------------------------------------------------------------------
// Constructors
Complex::Complex()
{
  real      = 0;
  imaginary = 0;
}
Complex::Complex( double new_real )
{
  real      = new_real;
  imaginary = 0;
}
Complex::Complex( double new_real, double new_imaginary )
{
  real      = new_real;
  imaginary = new_imaginary;
}
// Accessors
double Complex::getReal()
{
  return real;
}
double Complex::getIm()
{
  return imaginary;
}
// Setters
void Complex::setReal( double new_real )
{
  real = new_real;
  return;
}
void Complex::setIm( double new_imaginary )
{
  imaginary = new_imaginary;
  return;
}
void Complex::setNumber( double new_real, double new_imaginary )
{
  setReal( new_real );
  setIm( new_real );
  return;
}
// Conjugate
Complex Complex::getConjugate() 
{
  return Complex( real, -imaginary );
}
Complex* getConjugate( Complex* complex_vector, int size )
{
  Complex* conjugate_vector = (Complex*) malloc( size*sizeof(Complex) );
  for( int i=0; i<size; i++ )
  {
    conjugate_vector[i] = complex_vector[i].getConjugate();
  }
  return conjugate_vector;
}
// Basic Operations
Complex sum( Complex complex1, Complex complex2 )
{
  return Complex( complex1.getReal() + complex2.getReal(), complex1.getIm() + complex2.getIm() );
}
Complex multiply( Complex complex1, Complex complex2 )
{
  return Complex( complex1.getReal()*complex2.getReal() - complex1.getIm()*complex2.getIm(),      \
                  complex1.getReal()*complex2.getIm()   + complex1.getIm()*complex2.getReal() );
}
//----------------------------------------------------------------------

// Complex Exponential form
//-----------------------------------------------------------------
// Constructors
ComplexExp::ComplexExp()
{
  mods = 0;
  phase   = 0;
}
ComplexExp::ComplexExp( double new_modulus )
{
  mods    = new_modulus;
  phase   = 0;
}
ComplexExp::ComplexExp( double new_modulus, double new_phase )
{
  mods    = new_modulus;
  phase   = new_phase;
}
// Accessors
double ComplexExp::getModulus()
{
  return mods;
}
double ComplexExp::getPhase()
{
  return phase;
}
// Setters
void ComplexExp::setModulus( double new_mod )
{
  mods = new_mod;
}
void ComplexExp::setPhase( double new_phase )
{
  phase = new_phase;
}
void ComplexExp::setNumber( double new_mod, double new_phase )
{
  mods  = new_mod;
  phase = new_phase;
}
//-----------------------------------------------------------------