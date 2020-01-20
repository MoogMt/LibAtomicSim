#include "atom.h"

//==========
// Position
//======================================================================
std::vector<double> getPosition( AtomList & atoms , int index )
{
  std::vector<double> position;
  position.push_back( atoms.x[ index ] );
  position.push_back( atoms.y[ index ] );
  position.push_back( atoms.z[ index ] );
  return position;
}
//======================================================================
//===========
// DISTANCES
//======================================================================
double distanceAtoms(Atom i, Atom j)
// Returns the distance between two atoms
{
  double x2 = i.x - j.x;
  double y2 = i.y - j.y;
  double z2 = i.z - j.z;
  double dist=sqrt( pow(x2,2.0) + pow(y2,2.0) + pow(z2,2.0) );
  return dist;
}
//-------------------------------------------------------------------------
std::vector<double> distanceFromPoint( AtomList atom_list , std::vector<double> point )
{
  std::vector<double> r; r.assign( atom_list.x.size(), 0);
  std::vector<double> x2; x2 = square( distance( atom_list.x , point[0] ) );
  std::vector<double> y2; y2 = square( distance( atom_list.y , point[1] ) );
  std::vector<double> z2; z2 = square( distance( atom_list.z , point[2] ) );
  return squaroot( addVector( addVector( addVector( r , x2 ) , y2 ) , z2 ) );
}
//=====================================================================


//=======
// MOVE
//===========================================================================================
std::vector<Atom> compressAtoms( std::vector<Atom> atoms, double frac_a , double frac_b , double frac_c )
{
  for ( int i=0 ; i < atoms.size() ; i++ )
    {
      atoms[i].x *= frac_a;
      atoms[i].y *= frac_b;
      atoms[i].z *= frac_c;
    }
  return atoms;
}
//==========================================================================================

//=======
// IO
//======================================================================================
//-------
// PRINT
//--------------------------------------------------------------------------------------
void printAtoms( std::vector<Atom> atoms )
// Prints all informations on a given atoms in the console
{
  for( int i=0; i < atoms.size(); i++ )
    {
      std::cout << "ATOM " << i+1 << " | " << atoms[i].name << " | "  << atoms[i].x << " " << " " << atoms[i].y << " " << atoms[i].z << std::endl;
    }
  return ;
}
//--------
// WRITE
//--------------------------------------------------------------------------------------
void writePositions( std::ofstream & file , std::vector<Atom> atoms, std::string specie )
{
  for ( int i=0 ; i < atoms.size() ; i++ )
    {
      if ( atoms[i].name == specie )
	{
	  file << atoms[i].x << " ";
	  file << atoms[i].y << " ";
	  file << atoms[i].z << " ";
	  file << std::endl;
	}
    }
  file << std::endl;
  return;
}
//--------------------------------------------------------------------------------------
void writePositions( std::ofstream & file , AtomList atom_list , std::string specie )
{
  for ( int i=0 ; i < atom_list.x.size() ; i++ )
    {
      if ( atom_list.names[i] == specie )
	{
	  file << atom_list.x[i] << " ";
	  file << atom_list.y[i] << " ";
	  file << atom_list.z[i] << " ";
	  file << std::endl;
	}
    }
  file << std::endl;
  return;
}
//======================================================================================

