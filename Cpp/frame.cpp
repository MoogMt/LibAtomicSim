#include <assert.h> 
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>

#include "utils.h"
#include "frame.h"

using namespace std;

// Constructors      
Frame::Frame()
{
  number      = 0;
  names       = initWordVector(0);
  positions   = zeroMatrix(3,number);
  cell_matrix = identityMatrix(3);
}
Frame::Frame( int new_number )
{
  number      = new_number;
  names       = initWordVector( new_number ) ;
  positions   = zeroMatrix( 3, new_number );
  cell_matrix = initCellMatrix();
}
Frame::Frame( int new_number, std::string* new_names, double* new_positions )
{
  number    = new_number;
  names     = initWordVector( new_number );
  positions = zeroMatrix(3,new_number);
  cell_matrix=initCellMatrix();
  
  for ( int i=0; i < new_number; i++ )
  {
    int offset=3*i;
    names[i]=new_names[i];
    for (int j=0; j<3 ; j++ )
    {
      positions[offset+j]=new_positions[offset+j];
    }
  }
}

// Getters
// -> Number of atoms
int Frame::getNumberAtoms()
{
  return number; 
}
// -> Names of atoms
std::string* Frame::getNames()
{
  return names;
}
// -> Positions of atoms
double* Frame::getPositions()
{
 return positions;
}
// -> Cell matrix of the frame
double* Frame::getCellMatrix()
{
  return cell_matrix;
}

// Setters
// -> Number of atoms
void Frame::setNumberAtoms( int new_number )
{
  number = new_number;
}
// -> Names
void Frame::setNames( std::string* new_names ) 
{
  // Loop over atoms
  for ( int i=0; i < this->number; i++ )
  {
    // Setting name of atom i
    names[i] = new_names[i];
  }
  return;
}
void Frame::setName( std::string new_name, int index_atom )
{
  try
  {
    if ( index_atom < this->number )
    {
      names[index_atom] = new_name; 
    }
    else
    {
      throw(index_atom);
    }
  }
  catch( int index_atom )
  {
    std::cout << "Error!" << std::endl;
    std::cout << "Atom index ot of bond: " << index_atom << std::endl;
    std::cout << "Max index: " << this->number << std::endl;
  }
}
void Frame::setNames( double* new_names, int* index_atoms, int number_new ) 
{
  // Loop over atoms
  for( int i=0; i<number_new; i++)
  {
    // Catch error if the index is larger than the array size
    try
    {
      if ( index_atoms[i] < this->number )
      {
        positions[index_atoms[i]] = new_names[i];
      }
      // Send error message with the index of atom as message
      else
      {
        throw(index_atoms[i]);
      }
    }
    // IF the error is caught, then send a message in prompt
    catch( int index)
    {
      std::cout << "Error!" << std::endl;
      std::cout << "Index out of bond: " << index << std::endl;
      std::cout << "Index max: " << this->number << std::endl;
    }
  }
  return;
}
// -> Position
void Frame::initPosition( double* new_positions, int new_number )
{
  positions=zeroMatrix(3,new_number);
  for( int i=0; i<new_number; i++)
  {
    int offset=3*i;
    for( int j=0; j<3; j++ )
    {
      try
      {
        if(offset+j)
        {
          positions[offset+j]=new_positions[offset+j];
        }
        else
        {
          throw( offset+j );
        }
      }
      catch( int index )
      {
        indexOutOfBondError( index, this->number );
        return;
      }
    }
  }
  return;
}
void Frame::setPosition( double* new_positions, int index_atom )
{
  try
  { 
    if( index_atom < this->number )
    {
      int offset = 3*index_atom;
      for( int j=1; j<3 ; j++ )
      {
        positions[ offset+j ] = new_positions[j];
      }
    }
    else
    {
      throw( index_atom );
    }
  }
  catch( int index )
  {
    indexOutOfBondError( index, this->number );
  }

  return;
}
void Frame::setPositions( double* new_positions, int* index_atoms, int number_new ) 
{
  // Loop over atoms
  for( int i=0; i<number_new; i++)
  {
    // Catch error if the index is larger than the array size
    try
    {
      if ( index_atoms[i] < this->number )
      {
        // Bookkeeping
        int offset  = 3*index_atoms[i];
        int offset2 = 3*i;
        // Loop over dimensions
        for( int j=0; j<3; j++ )
        {
          // Change position
          if( offset2+j < 3*number_new )
          {
            positions[offset+j]=new_positions[offset2+j];
          }
          else
          {
            throw( offset2+j );
          }
        }
      }
      // Send error message with the index of atom as message
      else
      {
        throw(index_atoms[i]);
      }
    }
    catch(int index)
    {
      indexOutOfBondError( index, this->number );
      return;
    }
  }
  return;
}
void Frame::setPositions( double* new_positions )
{
  for ( int i=0; i < 3*number ; i++ )
  {
    int offset = 3*i;
    for( int j=0; i < 3; j++ )
    {
      try
      {
        if( offset+j < number )
        {
          positions[offset+j]=new_positions[offset+j];
        }
        else
        {
          throw( offset+j);
        }
      }
      catch( int index )
      {
         indexOutOfBondError( index, this->number );
         return;
      }
    }
  }
  return;
}

// Modifyers
// addAtoms: add a list of atoms to the original Frame
void Frame::addAtoms( int number_to_add, std::string* new_names, double* new_positions )
{
  // Internal variables
  int offset = 0;    // Avoid recomputing index
  int offset2 = 0;   // Avoid recomputing index (2)
  int old_number=number;  // Store the original number of atoms in the Frame
  std::string old_names[old_number]; // Stores the names of the original atoms in the Frame
  double old_positions[3*old_number]; // Stores the positions of the original atoms in the Frame 

  // Update the number of atoms in the Frame
  number = old_number + number_to_add; 

  // Storing original atoms names and positions
  for ( int i=0; i<old_number; i++ )
  {
    // Bookkeeping index of atom i
    offset = 3*i;
    // Stores name of atom i
    old_names[i] = names[i];
    // Loop overdimension
    for ( int j=0; j<3; j++ )
    {
      // Stores position j of atom i 
      old_positions[offset+j] = positions[offset+j];
    }
  }

  // Reallocating the Frame to make room for the new atoms names and positions
  names=(std::string*) malloc( (number)*sizeof(std::string) ) ;
  positions = (double*) malloc( (3*number)*sizeof(double) );

  // Moving back the original atoms names and positions to the new Frame
  for ( int i=0; i<old_number; i++ )
  {
    // Bookkeeping index of atom i
    offset = 3*i;
    // Moves name for atom i 
    names[i] = old_names[i];
    // Loop over dimensions
    for ( int j=0; j<3; j++ )
    {
      // Moving positions j of atom i
      positions[offset+j] = old_positions[offset+j];
    }
  }

  // Adding the new atoms names and positions to the Frame
  for ( int i=0; i<number_to_add; i++ )
  {
    // Bookkeeping the index for a given atom
    offset  = 3*(i + old_number);
    offset2 = 3*i;
    // Moving name for atom i
    names[ i + old_number ] = new_names[i];
    // Loop over dimensions
    for ( int j=0; j<3; j++ )
    {
      // Moving positions j for atom i 
      positions[ offset + j ] = new_positions[ offset2+j ];
    }
  }

  return;
}

// Checkers
// Print Atom Informations
void Frame::printAtoms() 
{
  for ( int i=0; i<(*this).getNumberAtoms(); i++ )
  {
    std::cout << (*this).getNames()[i] << " ";
    for ( int j=0; j<3; j++ )
    {
      std::cout << (*this).getPositions()[ (i*3) + j ] << " ";
    }
    std::cout << std::endl;
  }
  return;
}


