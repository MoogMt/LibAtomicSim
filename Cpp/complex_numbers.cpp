#include <assert.h> 
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>

#include "complex_numbers.h"

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

// Get accessors
double Complex::getReal()
{
  return real;
}
double Complex::getIm()
{
  return imaginary;
}

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

/*
Complex operator+( Complex complex1, Complex complex2 )
{
  Complex sum;

  return sum;
}*/


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

double ComplexExp::getModulus()
{
  return mods;
}
double ComplexExp::getPhase()
{
  return phase;
}

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