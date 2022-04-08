#include <assert.h> 
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

#include "general_io.h"

using namespace std;

bool fileExistsIO( ifstream &file_handler )
{
  if( file_handler.is_open() )
  {
    return true;
  }
  else
  {
    return false;
  }  
}
bool fileExists( std::string file_path )
{
  bool status=false;
  std::ifstream file_handler(file_path);
  status=fileExistsIO(file_handler);
  file_handler.close();
  return status;
}                                                                                                                       

int getNumberLineIO( ifstream &file_handler )
{
  int number_line=0;
  std::string str;
  while( getline(file_handler, str) )
  {
    number_line++;
  }

  if( file_handler.eof() )
    return number_line; 
  else
  {
    return -1;
  }
}
int getNumberLine( std::string file_path )
{
  int number_line=0;
  std::ifstream file_handler(file_path);
  number_line = getNumberLineIO(file_handler);
  file_handler.close();
  return number_line;
}

int countElement( std::string line, char delim )
{
  int counter=0;
  std::string buffer;
  stringstream line_stream(line);
  while( getline( line_stream, buffer, delim ) )
  {
    counter++;
  }
  return counter;
}

std::string* splitString( std::string target_string, char delim )
{
  int nb_element = countElement( target_string, delim );
  std::string string_element;
  stringstream string_stream(target_string);
  std::string split_string[nb_element];
  getline( string_stream, string_element, delim );
  split_string[0] = string_element;
  return split_string;
}