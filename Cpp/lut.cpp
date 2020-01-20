#include "lut.h"

//=======
// TESTS
//======================================================================
bool testType ( const AllTypeLUT lut_all, const int atom_index ,  const std::string specie )
{
  if ( lut_all.type_name[ atom_index ] == specie ) return true;
  else return false;
}
//---------------------------------------------------------------
bool testType ( const AllTypeLUT lut_all, const int atom_index ,  const int type_index )
{
  if ( lut_all.type_index[ atom_index ] == type_index )  return true;
  else return false;
}
//---------------------------------------------------------------
bool testType( const TypeLUT lut , const std::string specie )
{
  if ( lut.name == specie ) return true;
  else return false;
}
//---------------------------------------------------------------
bool testType( const TypeLUT lut , const int index )
{
  if ( lut.index == index ) return true;
  else return false;
}
//---------------------------------------------------------------
bool typeExists( const std::vector<TypeLUT> list , const std::string specie )
{
  for ( int i = 0 ; i < list.size() ; i++ )
    {
      if ( testType( list[i], specie ) ) return true;
    }
  return false;
}
//---------------------------------------------------------------
bool typeExists( const std::vector<TypeLUT> list , const int index )
{
  for ( int i = 0 ; i < list.size() ; i++ )
    {
      if ( testType( list[i], index ) ) return true;
    }
  return false;
}
//---------------------------------------------------------------
bool typeExists( const AllTypeLUT lut , const std::string specie )
{
  if ( typeExists( lut.types, specie ) ) return true;
  else return false;
}
//---------------------------------------------------------------
bool typeExists( const AllTypeLUT lut , const int type_index )
{
  if ( typeExists( lut.types, type_index ) ) return true;
  else return false;
}
//======================================================================

//======================================================================
bool checkNames( std::vector<std::string> names )
{
  for ( int i=0 ; i < names.size()-1 ; i++ )
    {
      for ( int j=i+1 ;j < names.size() ; j++ )
	{
	  if ( names[i] == names[j] ) return false;
	}
    }
  return true;
}
//======================================================================

//==================
// ADD ATOM TO LUT
//======================================================================
void addAtom2LUT( std::vector<TypeLUT>  & list , const std::string name , const int index  )
{
  for ( int i=0 ; i < list.size() ; i++ )
    {
      if ( name == list[i].name )
	{
	  list[i].atom_index.push_back( index );
	  return;
	}
    }
  list.push_back( { name, list.size() , initVector( index ) } );
  return;
}
//--------------------------------------------------------------
void addAtom2LUT( std::vector<TypeLUT> & list , const Atom atom )
{
  for ( int i=0 ; i < list.size() ; i++ )
    {
      if ( atom.name == list[i].name )
	{
	  list[i].atom_index.push_back( atom.index );
	  return;
	}
    }
  list.push_back( { atom.name , list.size() , initVector( atom.index ) } );
  return ;
}
//--------------------------------------------------------------
void addAtom2LUT( AllTypeLUT & list , const std::string name , const int index  )
{
  list.type_name.push_back( name );
  list.type_index.push_back( index );
  addAtom2LUT( list.types , name , index );
}
//--------------------------------------------------------------
void addAtom2LUT( AllTypeLUT & list , const Atom atom )
{
  list.type_name.push_back( atom.name );
  bool check = false;
  for ( int i=0 ; i < list.types.size() ; i++ )
    {
       if ( list.types[i].name == atom.name )
	{
	  list.type_index.push_back( list.types[i].index );
	  check=true;
	  break;
	}
    }
  if ( ! check )
    {
      list.type_index.push_back( list.types.size()  );
    }
  addAtom2LUT( list.types , atom.name , atom.index );
}
//======================================================================

//==========
// MAKE LUT
//======================================================================
TypeLUT makeLUT( const std::string name , const int index )
{
  std::vector<int> atom_index;
  return { name, index, atom_index };
}
//--------------------------------------------------------------
TypeLUT makeLUT( const std::string name , const int index , const std::vector<int> atom_index )
{
  return { name , index , atom_index };
}
//--------------------------------------------------------------------------------------
void makeLUT( std::vector<TypeLUT> & lut_list , const std::vector<std::string> names )
{
  std::vector<int> atom_index;
  for ( int i=0 ; i < names.size() ; i++ )
    {
      lut_list.push_back( makeLUT( names[i] , i ) );
    }
  return ;
}
//-------------------------------------------------------------------------------------
std::vector<TypeLUT> makeLUT( const std::vector<Atom> atoms )
{
  std::vector<TypeLUT> list_lut;
  for ( int i=0 ; i < atoms.size() ; i++ )
    {
      addAtom2LUT( list_lut, atoms[i] );
    }
  return list_lut;
}
//--------------------------------------------------------------------------------------
void makeLUT( AllTypeLUT & lut_list , const std::vector<std::string> names )
{
  makeLUT( lut_list.types , names);
  return;
}
//=======================================================================================

//============
// GET SPECIE
//==================================================
std::vector<int> getSpecieIndex( const AllTypeLUT & lut_maj, std::string specie )
{
  std::vector<int> index;
  for ( int i=0 ; i < lut_maj.types.size(); i++ )
    {
      if ( lut_maj.types[i].name == specie )
	{
	  return lut_maj.types[i].atom_index;
	}
    }
  if ( index.size() == 0 )
    {
      std::cout << "Chemical specie " << specie << " was not found in the box." << std::endl;
    }
  return index;
}
//==================================================
