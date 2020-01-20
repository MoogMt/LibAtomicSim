#include "histogram.h"

//========
// BINS
//======================================================================
Bin emptyBin()
{
  Bin bin = { 0, 0, 0};

  return bin;
}
//-----------------------------------------------------------
Bin makeBinReal( double begin , double end , int value )
{
  Bin bin = { begin , end , value };
  return bin;
}
//======================================================================

//========
// BINS
//======================================================================
double center( const Bin bin )
{
  return (bin.end+bin.begin)/2.0;
}
//-----------------------------------------------------------
bool overlap( const Bin bin1, const Bin bin2 )
{
  if ( max(bin1.begin,bin2.begin) <  min(bin1.end,bin2.end) )
    {
      return true;
    }
  else
    {
      return false;
    }
}
//-----------------------------------------------------------
void fillBin( Bin &bin , const std::vector<double> data )
{
  for ( int i=0 ; i < data.size() ; i++ )
    {
      if ( data[i] > bin.begin && data[i] < bin.end )
	{
	  bin.value++;
	}
    }
  return ;
}
//-----------------------------------------------------------
Bin makeBin( const double bin_min, const double bin_max, const std::vector<double> data)
{
  Bin bin = { bin_min , bin_max , 0 };
  fillBin( bin , data );
  return bin;
}
//-----------------------------------------------------------
Bin addBinsMin( const Bin bin1, const Bin bin2)
{
  Bin bin_final;
  if ( overlap(bin1,bin2) )
    {
      bin_final.begin = max(bin1.begin,bin2.begin);
      bin_final.end   = min(bin1.end,bin2.end);
      bin_final.value = bin1.value + bin2.value;
      return bin_final;
    }
  else
    {
      return emptyBin();
    }
}
//-------------------------------------------------------------
Bin addBinsMax( const Bin bin1, const Bin bin2 )
{
  Bin bin_final = { min(bin1.begin,bin2.begin) , max(bin1.end,bin2.end) , bin1.value + bin2.value };
  return bin_final;
}
//===============================================================================================

//================
// MAKE HISTOGRAM
//===============================================================================================
std::vector<Bin> makeRegularHistogram( const std::vector<double> data , const double x_min , const double x_max , const int number_bins )
{
  std::vector<Bin> bins_hist;
  double delta = (x_max-x_min)/(double)(number_bins);
  if ( delta < 0 )
    {
      std::cout << "Error: Min of the interval to max!" << std::endl;
      return bins_hist;
    }
  for ( int i=0 ; i < number_bins ; i++ )
    {
      bins_hist.push_back( makeBin( x_min + i*delta , x_min + (i+1)*delta , data ) );
    }
  return bins_hist;
}
//--------------------------------------------------------------------------------------------
std::vector<Bin> makeRegularHistogram( const std::vector<double> data_x , const std::vector<double> data_y , const int number_bins )
{
  std::vector<Bin> bins_hist;
  double x_min = min(data_x);
  double x_max = max(data_x);
  double delta = (x_max-x_min)/(double)(number_bins);
  for ( int i=0 ; i < number_bins ; i++ )
    {
      bins_hist.push_back( makeBin( x_min + i*delta , x_min + (i+1)*delta , data_y ) );
    }
  return bins_hist;
}
//--------------------------------------------------------------------------------------------
std::vector<Bin> makeHistograms( const std::vector<double> data, const std::vector<double> bins_limits )
{
  std::vector<Bin> bins_hist;
  for( int i=0 ; i < bins_limits.size()-1 ; i++ )
    {
      bins_hist.push_back( makeBin( bins_limits[i] , bins_limits[i+1] , data ) );
    }
  return bins_hist;
}
//--------------------------------------------------------------------------------------------
std::vector<Bin> makeHistograms( const std::vector<double> data, std::vector<Bin> bins )
{
  for( int i=0 ; i < bins.size() ; i++ )
    {
      fillBin( bins[i] , data ) ;
    }
  return bins;
}
//===============================================================================================


//===============
// Modification
//===============================================================================================
std::vector<Bin> addHistograms( const std::vector<Bin> hist1 , const std::vector<Bin> hist2 )
{
  std::vector<Bin> sum_hist;
  for ( int i=0 ; i < hist1.size() ; i++ )
    {
      for ( int j=0 ; j  < hist2.size() ; j++ )
	{
	  if ( overlap( hist1[i] , hist2[j] ) )
	    {
	      sum_hist.push_back( addBinsMin( hist1[i] , hist2[j] ) );
	    }
	}
    }
  return sum_hist;
}
//===============================================================================================

//==============
// Total Value
//===============================================================================================
double getTotalValue( const std::vector<Bin> hist )
{
  double total_value = 0;
  for ( int i=0 ; i < hist.size() ; i++ )
    {
      total_value += hist[i].value;
    }
  return total_value;
}
//-----
// IO
//---------------------------------------------------------------------------------------
void writeHistogram( std::ofstream & file , const std::vector<Bin> hist )
{
  for ( int i = 0 ; i < hist.size() ; i++ )
    {
      file  << center(hist[i]) << " " << hist[i].value << std::endl;
    }
  return;
}
void writeHistogram( std::string file_name , const std::vector<Bin> hist )
{
  std::ofstream file ( file_name.c_str() ,  std::ios::out | std::ios::app );
  for ( int i=0 ; i < hist.size() ; i++ )
    {
      file << center(hist[i]) << " " << hist[i].value << std::endl;
    }
  file.close();
  return;
}
//==========================================================================================

//==========
// Bin Real
//==========================================================================================
BinReal emptyBinReal( )
{
  BinReal bin = { 0, 0, 0};
  return bin;
}
//-------------------------------------------------------------------------------------------------
BinReal makeBinReal( const double begin , const double end , const double value )
{
  BinReal bin = { begin , end , value };
  return bin;
}
//-------------------------------------------------------------------------------------------------
double center( const BinReal bin )
{
  return ( bin.begin + bin.end )/2;
}
//-------------------------------------------------------------------------------------------------
bool overlap( const BinReal bin1 , const BinReal bin2 )
{ 
  if ( max(bin1.begin,bin2.begin) <  min(bin1.end,bin2.end) )
    {
      return true;
    }
  else
    {
      return false;
    }
}
//-------------------------------------------------------------------------------------------------
void fillBin( BinReal &bin , const std::vector<double> data )
{
  for ( int i=0 ; i < data.size() ; i++ )
    {
      if ( data[i] > bin.begin && data[i] < bin.end )
	{
	  bin.value++;
	}
    }
  return ;
}
//-------------------------------------------------------------------------------------------------
BinReal makeBinReal( const double bin_min, const double bin_max, const std::vector<double> data)
{
  BinReal bin = { bin_min , bin_max , 0 };
  fillBin( bin , data );
  return bin;
}
//-------------------------------------------------------------------------------------------------
BinReal addBinsMin( const BinReal bin1, const BinReal bin2)
{
  BinReal bin_final;
  if ( overlap(bin1,bin2) )
    {
      bin_final.begin = max(bin1.begin,bin2.begin);
      bin_final.end   = min(bin1.end,bin2.end);
      bin_final.value = bin1.value + bin2.value;
      return bin_final;
    }
  else
    {
      return emptyBinReal();
    }
}
//-------------------------------------------------------------------------------------------------
BinReal addBinsMax( const BinReal bin1, const BinReal bin2 )
{
  BinReal bin_final = { min(bin1.begin,bin2.begin) , max(bin1.end,bin2.end) , bin1.value + bin2.value };
  return bin_final;
}
//================================================================================================


//======================
// Bin Real Histogram
//================================================================================================
double getTotalValue( const std::vector<BinReal> hist )
{
  double total_value = 0;
  for ( int i=0 ; i < hist.size() ; i++ )
    {
      total_value += hist[i].value;
    }
  return total_value;
}
//-------------------------------------------------------------------------------------------------
std::vector<BinReal> normalizeHistogram( const std::vector<Bin> hist )
{
  std::vector<BinReal> hist_real;
  double norm = getTotalValue( hist );
  for ( int i=0 ; i < hist.size() ; i++ )
    {
      hist_real.push_back( makeBinReal( hist[i].begin , hist[i].end , hist[i].value/norm ) );
    }
  return hist_real;
}
//-------------------------------------------------------------------------------------------------
std::vector<BinReal> normalizeHistogram( std::vector<BinReal> hist )
{
  double norm = getTotalValue( hist );
  for ( int i=0 ; i < hist.size() ; i++ )
    {
      hist[i].value = hist[i].value/norm;
    }
  return hist;
}
//-------------------------------------------------------------------------------------------------
void writeHistogram( std::ofstream & file , const std::vector<BinReal> hist )
{
  for ( int i=0 ; i < hist.size() ; i++ )
    {
      file << center(hist[i]) << " " << hist[i].value << std::endl;
    }
  return;
}
//-------------------------------------------------------------------------------------------------
void writeHistogram( std::string file_name , const std::vector<BinReal> hist )
{
  std::ofstream file ( file_name.c_str() ,  std::ios::out | std::ios::app );
  for ( int i=0 ; i < hist.size() ; i++ )
    {
      file << center(hist[i]) << " " << hist[i].value << std::endl;
    }
  file.close();
  return;
}
//-------------------------------------------------------------------------------------------------
bool checkSizeHists( const std::vector<BinReal> hist , const std::vector<BinReal> hist2 )
{
  if ( hist.size() == hist2.size() )
    {
      return true;
    }
  else
    {
      return false;
    }
}
//-------------------------------------------------------------------------------------------------
bool checkSizeHists( const std::vector< std::vector<BinReal> > hist_list )
{
  for ( int i=0 ; i < hist_list.size()-1 ; i++ )
    {
      for ( int j=0 ; j < hist_list.size() ; j++ )
	{
	  if ( ! checkSizeHists( hist_list[i] , hist_list[j] ) )
	    {
	      return false;
	    }
	}
    }
  return true;
}
//-------------------------------------------------------------------------------------------------
void writeBinReal( std::ofstream & file , const BinReal bin, const bool wcenter)
{
  if ( wcenter )
    {
      file << center(bin) <<  " " << bin.value << " ";
    }
  else
    {
      file << bin.value << " ";
    }
  return ;
}
//-------------------------------------------------------------------------------------------------
void writeBinRealCenter( std::ofstream & file, const BinReal bin )
{
  file << center(bin) << " ";
  return;
}
//-------------------------------------------------------------------------------------------------
void writeHistBinCenter( std::ofstream & file , const std::vector<BinReal> hist , const int index )
{
  writeBinRealCenter( file, hist[index] );
}  
//-------------------------------------------------------------------------------------------------
void writeHistBin( std::ofstream & file , const std::vector<BinReal> hist, const int index , const bool wcenter )
{
  writeBinReal( file , hist[index] , wcenter );
}
//-------------------------------------------------------------------------------------------------
void writeHistograms( std::ofstream & file , const std::vector< std::vector<BinReal> > hist_list )
{
  if ( checkSizeHists( hist_list ) )
    {
      std::cout << "Error: Some histograms do not have the same size." << std::endl;
      return;
    }
  for ( int i=0 ; i < hist_list[0].size() ; i++ )
    {
      writeHistBinCenter( file, hist_list[0] , i );
      for ( int j=0 ; j < hist_list.size() ; j++ )
	{
	  writeHistBin( file , hist_list[j] , i , false );
	}
      file << std::endl;
    }  
  return ;
}
//=========================================================================================

//=======================
// Integrating Histogram
//=========================================================================================
double average( const std::vector<BinReal> & histogram )
{
  double value =0;
  for ( int i=0; i < histogram.size() ; i++ )
    {
      value += center(histogram[i])*histogram[i].value;
    }
  return value;
}
//------------------------------------------------------------
double variance( const std::vector<BinReal> & histogram )
{
  double value  = 0;
  double value2 = 0;
  for ( int i=0; i < histogram.size() ; i++ )
    {
      value += center(histogram[i])*histogram[i].value;
      value2 += center(histogram[i])*center(histogram[i])*histogram[i].value;
    }
  return sqrt(value2-value*value);
}
//------------------------------------------------------------
double integrateHistogram( const std::vector<BinReal> & histogram )
{
  double value = 0;
  for ( int i=0 ; i < histogram.size() ; i++ )
    {
      value += histogram[i].value;
    }
  return value;
}
//------------------------------------------------------------
double integrateHistogram( const std::vector<BinReal> & histogram , const double start , const double end )
{
  double value = 0;
  for ( int i=0 ; i < histogram.size() ; i++ )
    {
      if ( histogram[i].begin > start && histogram[i].end < end )
	{
	  value += histogram[i].value;
	}
    }
  return value;
}
//=========================================================================================

//=====
// IO
//=========================================================================================
void readRegularHistogram( const std::string file_name , std::vector<Bin> & histogram )
{

  // Reading variables
  std::ifstream file( file_name.c_str() );
  std::string line, line1, line2;
  double step = 0;

  // First two step, to compute the step
  if( std::getline( file , line1 ) )
    {
      double step1 = 0, value1 = 0;
      std::istringstream it_string1(line1);
      if ( it_string1 >> step1 >> value1 ) 
	{
	  if( std::getline( file , line2 ) )
	    {
	      std::istringstream it_string2( line2 );
	      double step2 = 0 , value2 = 0;
	      if ( it_string2 >> step2 >> value2 )
		{
		  step = (step2 - step1)*0.5;
		  Bin bin1 = { step1 - step , step1 + step , value1 };
		  Bin bin2 = { step2 - step , step2 + step , value2 };
		  histogram.push_back( bin1 );
		  histogram.push_back( bin2 );
		}
	      else return;
	    }
	  else return;
	}
      else return;
    }
  else return;

  while( std::getline( file , line ) )
    {
      double step_loc, value;
      std::istringstream it_string( line );
      if ( it_string >> step_loc >> value )
	{
	  Bin bin = { step_loc - step , step_loc + step , value };
	  histogram.push_back( bin );
	}
      return;
    }
  return;
}
//-----------------------------------------------------------------------------------
void readRegularHistogram( const std::string file_name , std::vector<BinReal> & histogram )
{

  // Reading variables
  std::ifstream file( file_name.c_str() );
  std::string line, line1, line2;
  double step = 0;

  // First two step, to compute the step
  if( std::getline( file , line1 ) )
    {
      double step1 = 0, value1 = 0;
      std::istringstream it_string1(line1);
      if ( it_string1 >> step1 >> value1 ) 
	{
	  if( std::getline( file , line2 ) )
	    {
	      std::istringstream it_string2( line2 );
	      double step2 = 0 , value2 = 0;
	      if ( it_string2 >> step2 >> value2 )
		{
		  step = (step2 - step1)*0.5;
		  BinReal bin1 = { step1 - step , step1 + step , value1 };
		  BinReal bin2 = { step2 - step , step2 + step , value2 };
		  histogram.push_back( bin1 );
		  histogram.push_back( bin2 );
		}
	      else return;
	    }
	  else return;
	}
      else return;
    }
  else return;

  while( std::getline( file , line ) )
    {
      double step_loc, value;
      std::istringstream it_string( line );
      if ( it_string >> step_loc >> value )
	{
	  BinReal bin = { step_loc - step , step_loc + step , value };
	  histogram.push_back( bin );
	}
      return;
    }
  return;
}
//-----------------------------------------------------------------------------------
std::vector<Bin> readRegularHistogram( const std::string file_name )
{
  //-----------
  // Histogram
  //--------------------------------------------------------------
  std::vector<Bin> histogram;
  //--------------------------------------------------------------

  //-------------------
  // Reading variables
  //--------------------------------------------------------------
  std::ifstream file( file_name.c_str() );
  std::string line, line1, line2;
  double step = 0;
  //--------------------------------------------------------------

  //--------------------------------------------------------------
  // First two step, to compute the step
  //--------------------------------------------------------------
  if( std::getline( file , line1 ) )
    {
      double step1 = 0, value1 = 0;
      std::istringstream it_string1(line1);
      if ( it_string1 >> step1 >> value1 ) 
	{
	  if( std::getline( file , line2 ) )
	    {
	      std::istringstream it_string2( line2 );
	      double step2 = 0 , value2 = 0;
	      if ( it_string2 >> step2 >> value2 )
		{
		  step = (step2 - step1)*0.5;
		  Bin bin1 = { step1 - step , step1 + step , value1 };
		  Bin bin2 = { step2 - step , step2 + step , value2 };
		  histogram.push_back( bin1 );
		  histogram.push_back( bin2 );
		}
	      else return histogram;
	    }
	  else return histogram;
	}
      else return histogram;
    }
  else return histogram;
  //--------------------------------------------------------------

  //-------------------------------
  // Reading rest of the histogram
  //--------------------------------------------------------------
  while( std::getline( file , line ) )
    {
      double step_loc, value;
      std::istringstream it_string( line );
      if ( it_string >> step_loc >> value )
	{
	  Bin bin = { step_loc - step , step_loc + step , value };
	  histogram.push_back( bin );
	}
      else return histogram;
    }
  //--------------------------------------------------------------
  
  // Return histogram
  return histogram;
}
//-----------------------------------------------------------------------------------
std::vector<BinReal> readRegularHistogramReal( const std::string file_name )
{
  //-----------
  // Histogram
  //--------------------------------------------------------------
  std::vector<BinReal> histogram;
  //--------------------------------------------------------------

  //-------------------
  // Reading variables
  //--------------------------------------------------------------
  std::ifstream file( file_name.c_str() );
  std::string line, line1, line2;
  double step = 0;
  //--------------------------------------------------------------

  //--------------------------------------------------------------
  // First two step, to compute the step
  //--------------------------------------------------------------
  if( std::getline( file , line1 ) )
    {
      double step1 = 0, value1 = 0;
      std::istringstream it_string1(line1);
      if ( it_string1 >> step1 >> value1 ) 
	{
	  if( std::getline( file , line2 ) )
	    {
	      std::istringstream it_string2( line2 );
	      double step2 = 0 , value2 = 0;
	      if ( it_string2 >> step2 >> value2 )
		{
		  step = (step2 - step1)*0.5;
		  BinReal bin1 = { step1 - step , step1 + step , value1 };
		  BinReal bin2 = { step2 - step , step2 + step , value2 };
		  histogram.push_back( bin1 );
		  histogram.push_back( bin2 );
		}
	      else return histogram;
	    }
	  else return histogram;
	}
      else return histogram;
    }
  else return histogram;
  //--------------------------------------------------------------

  //-------------------------------
  // Reading rest of the histogram
  //--------------------------------------------------------------
  while( std::getline( file , line ) )
    {
      double step_loc, value;
      std::istringstream it_string( line );
      if ( it_string >> step_loc >> value )
	{
	  BinReal bin = { step_loc - step , step_loc + step , value };
	  histogram.push_back( bin );
	}
      else return histogram;
    }
  //--------------------------------------------------------------
  
  // Return histogram
  return histogram;
}
//=========================================================================================

