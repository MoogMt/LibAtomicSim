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
double Complex::re()
{
  return real;
}
double Complex::im()
{
  return imaginary;
}
// Setters
void Complex::re( double new_real )
{
  real = new_real;
  return;
}
void Complex::im( double new_imaginary )
{
  imaginary = new_imaginary;
  return;
}
void Complex::set( double new_real, double new_imaginary )
{
  re( new_real );
  im( new_real );
  return;
}
// Conjugate
Complex Complex::conj() 
{
  return Complex( real, -imaginary );
}
Complex* conj( Complex* complex_vector, int size )
{
  Complex* conjugate_vector = (Complex*) malloc( size*sizeof(Complex) );
  for( int i=0; i<size; i++ )
  {
    conjugate_vector[i] = complex_vector[i].conj();
  }
  return conjugate_vector;
}
Complex Complex::mod( Complex complex )
{
  return sqrt( real*real + imaginary*imaginary );
}
// Basic Operations
Complex sum( Complex complex1, Complex complex2 )
{
  return Complex( complex1.re() + complex2.re(), complex1.im() + complex2.im() );
}
Complex multiply( Complex complex1, Complex complex2 )
{
  return Complex( complex1.re()*complex2.re() - complex1.im()*complex2.im(),      \
                  complex1.re()*complex2.im() + complex1.im()*complex2.re() );
}
//----------------------------------------------------------------------

// Complex Exponential form
//-----------------------------------------------------------------
// Constructors
ComplexExp::ComplexExp()
{
  modulus = 0;
  phase   = 0; // in radians
}
ComplexExp::ComplexExp( double new_mod )
{
  modulus = new_mod;
  phase   = 0;
}
ComplexExp::ComplexExp( double new_mod, double new_phase )
{
  modulus = new_mod;
  phase   = new_phase;
}
// Accessors
double ComplexExp::mod()
{
  return modulus;
}
double ComplexExp::phi()
{
  return phase;
}
// Setters
void ComplexExp::mod( double new_mod )
{
  modulus = new_mod;
}
void ComplexExp::phi( double new_phase )
{
  phase = new_phase;
}
void ComplexExp::set( double new_mod, double new_phase )
{
  modulus = new_mod;
  phase   = new_phase;
}
// Get Real and Imaginary part
double ComplexExp::re()
{
  return modulus*cos(phase);
}
double ComplexExp::im()
{
  return modulus*sin(phase);
}
//-----------------------------------------------------------------