#include <assert.h> 
#include <iostream>
#include <ostream>
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
void Complex::set( Complex complex )
{
  re( complex.re() );
  im( complex.im() );
  return;
}
// Operators
// -> Conjugate
Complex Complex::getConj() 
{
  return Complex( real, -imaginary );
}
void Complex::setConj()
{
  this->im( -this->im() );
  return;
}
Complex* getConj( Complex* complex_vector, int size )
{
  Complex* conjugate_vector = (Complex*) malloc( size*sizeof(Complex) );
  for( int i=0; i<size; i++ )
  {
    conjugate_vector[i] = complex_vector[i].getConj();
  }
  return conjugate_vector;
}
void setConj( Complex* complex_vector, int size )
{
  for( int i=0; i<size; i++ )
  {
    complex_vector[i] = complex_vector[i].getConj();
  }
  return;
}
// -> Modulus
double Complex::mod()
{
  return sqrt( this->re()*this->re() + this->im()*this->im() );
}
// -> All arithmetic operations (except modulus)
Complex Complex::operator+( Complex complex2 )
{
  return Complex( this->re() - complex2.re(), this->im()-complex2.im() );
}
Complex Complex::operator+( double real )
{
  return Complex( this->re() + real , this->im() );
}
void Complex::operator+=( Complex complex )
{
  this->set( this->re() + complex.re(), this->im() + complex.im() );
  return;
}
void Complex::operator+=( double real )
{
  this->set( this->re() + real, this->im() );
  return;
}
Complex Complex::operator-( Complex complex2 )
{
  return Complex( this->re() - complex2.re(), this->im()-complex2.im() );
}
Complex Complex::operator-( double real )
{
  return Complex( this->re() - real, this->im() );
}
Complex operator-( Complex complex )
{
  return Complex( - complex.re(), -complex.im() );
}
void Complex::operator-=( Complex complex )
{
  this->set( this->re() - complex.re(), this->im() - complex.im() );
  return;
}
void Complex::operator-=( double real )
{
  this->set( this->re() - real, this->im() );
  return;
}
Complex Complex::operator*( Complex complex2 )
{
  return Complex( this->re()*complex2.re() - this->im()*complex2.im(), this->im()*complex2.re() + this->re()*complex2.im() );
}
Complex Complex::operator*( double real )
{
  return Complex( real*this->re(), real*this->im() );
}
void Complex::operator*=( Complex complex )
{
  this->set( this->re()*complex.re() - this->im()*complex.im(), this->im()*complex.re() + this->re()*complex.im() );
  return;
}
void Complex::operator*=( double real )
{
  this->set( real*this->re(), real*this->im() );
  return;
}
Complex Complex::operator/( Complex complex2 )
{
  double mod2squared = pow( complex2.mod(),2);
  return Complex( (this->re() + complex2.re())/mod2squared, (this->im() - complex2.im())/mod2squared );
}
Complex Complex::operator/( double real )
{
  try
  {
    if( abs(real) > 0.0000000000000000001 )
    {
      return Complex( this->re()/real, this->im()/real );
    }
    else
    {
      throw( real );
    }
  }
  catch( double value )
  {
    std::cout << "Can't divide a number by 0" << std::endl;
    std::cout << "Value: " << real << std::endl;
    return Complex(0,0);
  }
}
void Complex::operator/=( Complex complex )
{

  double smaller = 0.00000000000000001;
  try 
  {
    if ( abs( complex.re() > smaller ) && abs( complex.im() ) > smaller  )
    {
      double mod2squared = pow( complex.mod(),2);
      this->set( (this->re() + complex.re() )/mod2squared, ( this->im() - complex.im() )/mod2squared );
      return;
    }
    else
    {
      throw( smaller );
    }
  }
  catch( double value )
  {
    std::cout << "Can't divide a number by 0" << std::endl;
    std::cout << "Value: " << real << std::endl;
    this->set(0,0);
    return;
  }
}
void Complex::operator/=( double real )
{
  double smaller = 0.0000000000001;
  try
  {
    if( abs(real) > smaller )
    {
      this->set( this->re()/real, this->im()/real );
      return;
    }
    else
    {
      throw( real );
    }
  }
  catch( double value )
  {
    std::cout << "Can't divide a number by 0" << std::endl;
    std::cout << "Value: " << real << std::endl;
    this->set(0,0);
    return;
  }
}
// -> Output flux operator
std::ostream& operator << ( std::ostream& os, const Complex& complex )
{
  double smallest = 0.000000000000001 ; 
  if( abs(complex.real) > smallest )
  {
    os << complex.real;
    if( abs( complex.imaginary ) > smallest )
    {
      os << " + " << complex.imaginary << "i";
    }
  }
  else
  {
    os << complex.imaginary << "i";
  }
  return os;
}
//----------------------------------------------------------------------

// Complex Exponential form
//-----------------------------------------------------------------
// Constructors
ComplexExp::ComplexExp()
{
  rho = 0;
  phi = 0; // in radians
}
ComplexExp::ComplexExp( double new_mod )
{
  rho = new_mod;
  phi = 0;
}
ComplexExp::ComplexExp( double new_mod, double new_phase )
{
  rho = new_mod;
  phi = new_phase;
}
// Accessors
double ComplexExp::mod()
{
  return rho;
}
double ComplexExp::phase()
{
  return phi;
}
// Setters
void ComplexExp::mod( double new_mod )
{
  rho = new_mod;
}
void ComplexExp::phase( double new_phase )
{
  phi = new_phase;
}
void ComplexExp::set( double new_mod, double new_phase )
{
  rho = new_mod;
  phi = new_phase;
}
// Get Real and Imaginary part
double ComplexExp::re()
{
  return rho*cos(phi);
}
double ComplexExp::im()
{
  return rho*sin(phi);
}
// Operators
ComplexExp ComplexExp::operator*( ComplexExp complex )
{
  return ComplexExp( this->mod()*complex.mod(), fmod( this->phase()+complex.phase(), 2*M_PI ) );
}
ComplexExp ComplexExp::operator*( double real )
{
  return ComplexExp( this->mod()*real, this->phase() );
}
ComplexExp ComplexExp::operator/( ComplexExp complex )
{
  return ComplexExp( this->mod()/complex.mod(), fmod( this->phase()-complex.phase(), 2*M_PI ) );
}
ComplexExp ComplexExp::operator/( double real )
{
  return ComplexExp( this->mod()/real, this->phase() );
}
//-----------------------------------------------------------------