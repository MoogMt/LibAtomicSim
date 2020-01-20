#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>  
#include <stdlib.h>    
#include <math.h>
#include <sstream>

#include "utils.h"

//------
// BINS
//----------------------
struct Bin
{
  double begin, end; // Inner and Outer bound of the bin
  int value;         // Number of elements inside the bin
};
//----------------------
struct BinReal
{
  double begin, end; // Inner and Outer bound of the bin
  double value;         // Number of elements inside the bin
};
//----------------------

//===========
// FUNCTIONS
//==============================================================================

//------
// Bins
//------------------------------------------------------------------------------
Bin emptyBin();
double center( const Bin bin );
Bin makeBinReal( const double begin , const double end , const int value );
bool overlap( const Bin bin1, const Bin bin2 );
void fillBin( Bin &bin , const std::vector<double> data );
Bin makeBin( const int bin_min, const int bin_max, const std::vector<double> data );
Bin addBinsMin( const Bin bin1, const Bin bin2 );
Bin addBinsMax( const Bin bin1, const Bin bin2 );
//------------------------------------------------------------------------------

//----------------
// Histograms
//------------------------------------------------------------------------------
// Regular histogram
//------------------------------------------------------------------------------
std::vector<Bin> makeRegularHistogram( const std::vector<double> data , const double x_min , const double x_max , const int number_bins );
std::vector<Bin> makeRegularHistogram( const std::vector<double> data_x , const std::vector<double> data_y , const int number_bins );
// Custom histogram
//------------------------------------------------------------------------------
std::vector<Bin> makeHistograms( const std::vector<double> data, const std::vector<double> bins_limits );
std::vector<Bin> makeHistograms( const std::vector<double> data, const std::vector<Bin> bins );
std::vector<Bin> addHistograms( const std::vector<Bin> hist1 , const std::vector<Bin> hist2 );
double getTotalValue( const std::vector<Bin> hist );
void writeHistogram( std::ofstream & file , const std::vector<Bin> hist );
void writeHistogram( std::string file_name , const std::vector<Bin> hist );
//------------------------------------------------------------------------------

//-----------
// Bin Real
//------------------------------------------------------------------------------
BinReal emptyBinReal( );
double center( const BinReal bin );
bool overlap( const BinReal bin1 , const BinReal bin2 );
BinReal makeBinReal( const double begin , const double end , const double value );
void fillBin( BinReal &bin , const std::vector<double> data );
BinReal makeBinReal( const double bin_min, const double bin_max, const std::vector<double> dat);
BinReal addBinsMin( const BinReal bin1, const BinReal bin2 );
BinReal addBinsMax( const BinReal bin1, const BinReal bin2 );
//------------------------------------------------------------------------------

//-----------------
// Real Histogram
//------------------------------------------------------------------------------
double getTotalValue( const std::vector<BinReal> hist );
std::vector<BinReal> normalizeHistogram( const std::vector<Bin> hist );
std::vector<BinReal> normalizeHistogram( const std::vector<BinReal> hist );
void writeHistogram( std::ofstream & file , const std::vector<BinReal> hist );
void writeHistogram( const std::string file_name , const std::vector<BinReal> hist );
bool checkSizeHists( const std::vector<BinReal> hist , const std::vector<BinReal> hist2 );
bool checkSizeHists( const std::vector< std::vector<BinReal> > hist_list );
void writeBinReal( std::ofstream & file , const BinReal bin, const bool wcenter);
void writeBinRealCenter( std::ofstream & file, const BinReal bin );
void writeHistBinCenter( std::ofstream & file , const std::vector<BinReal> hist , const int index );
void writeHistBin( std::ofstream & file , const std::vector<BinReal> hist, const int index , const bool wcenter );
void writeHistograms( std::ofstream & file , const std::vector< std::vector<BinReal> > hist_list );
//------------------------------------------------------------------------------

//-----------------------
// Integrate Histograms
//-------------------------------------------------------------------------------------
double average( const std::vector<BinReal> & histogram );
double variance( const std::vector<BinReal> & histogram );
double integrateHistogram( const std::vector<BinReal> & histogram );
double integrateHistogram( const std::vector<BinReal> & histogram , const double start , const double end );
//-------------------------------------------------------------------------------------

//----
// IO
//-------------------------------------------------------------------------------------
void readRegularHistogram( const std::string file_name , std::vector<Bin> & histogram );
void readRegularHistogram( const std::string file_name , std::vector<BinReal> & histogram );
std::vector<Bin> readRegularHistogram( const std::string file_name );
std::vector<BinReal> readRegularHistogramReal( const std::string file_name );
//-------------------------------------------------------------------------------------

//======================================================================================

#endif 
