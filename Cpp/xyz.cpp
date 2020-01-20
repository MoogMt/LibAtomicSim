#include "xyz.h"

//=====================
// Get number of Steps
//=================================================================================
int getStepXYZ( std::ifstream & file )
{
  // Vars
  //--------------------------------------
  int nb_lines=0;    // Number of lines
  int nb_atoms=0;    // Number of atoms
  std::string line;  // Line
  //--------------------------------------
  
  // Get number of atoms
  //------------------------------------------
  if ( std::getline( file, line ) )
    {
      std::istringstream it_string(line);
      if ( ! ( it_string >> nb_atoms )) return -2;
    }
  else return -1;
  //------------------------------------------

  // Rewind
  file.seekg (0, file.beg);

  // Read whole file to count number of lines
  while ( std::getline( file, line ) ) nb_lines++;

  // Rewind
  file.seekg (0, file.beg);
  
  // Returns number of steps
  return nb_lines/(double)(nb_atoms);
}
//=================================================================================

//======
// READ
//=================================================================================
bool readStepXYZfast( std::ifstream & file , AtomList & atoms , AllTypeLUT & lut_list , const bool same_type , const bool verbose)
// Reads a step of XYZ file and build atom list and LUT for types
{
  //----------------------
  // Physical Parameters
  //---------------------------------------------------------------
  int nb_atoms = 0 ; // Numbers
  std::string line ; // Contains one line
  //---------------------------------------------------------------

  //------------------------------------------
  // Getting the number of atoms in the step
  //-----------------------------------------------------------------------------------
  if ( std::getline( file , line ) )
    {
      //std::cout << "line: " << line << std::endl;
      std::istringstream it_string(line);
      if ( ! ( it_string >> nb_atoms ) )
	{
	  if ( verbose ) std::cout << "Problem with file format at line 1." << std::endl;
	  return false;
	}
      //std::cout << "nb_atoms: " << nb_atoms << std::endl;
    }
  else return false;
  //-----------------------------------------------------------------------------------
  
  //-------------
  // Filling LUT
  //---------------------------------------------------------------
  bool fill_lut;
  if ( lut_list.type_name.size() == 0 ) fill_lut=true;
  else fill_lut = false;
  //---------------------------------------------------------------

  //-----------------------------
  // Initialization of variables
  //-----------------------------------------------------------------------------------
  if ( atoms.names.size() != nb_atoms )
    {
      atoms.names.assign( nb_atoms , "" );
      atoms.x.assign( nb_atoms, 0 );
      atoms.y.assign( nb_atoms, 0 );
      atoms.z.assign( nb_atoms, 0 );
      atoms.index.assign( nb_atoms , 0 );
    }
  //-----------------------------------------------------------------------------------

  //---------------
  // Detecting EOF
  //-----------------------------------------------------------------------------------
  if ( nb_atoms == 0 )
    {
      return false;
    }
  //-----------------------------------------------------------------------------------
  
  //-------------------------
  // Jumps next comment line
  //----------------------------------------------------------------------------------
  if( ! std::getline( file , line ) )
    {
      if ( verbose && line != "" )
	{
	  std::cout << "Problem with file format at line 2" << std::endl;
	}
      return false;
    }
  //----------------------------------------------------------------------------------

  //----------------------
  // Reads atoms by atoms
  //-------------------------------------------------------------------------------------
  int i = 0;
  while ( i < nb_atoms )
    {
      if ( std::getline( file , line ) )
	{
	  std::istringstream it_string(line);
	  if ( ! ( it_string >> atoms.names[i] >> atoms.x[i] >> atoms.y[i] >> atoms.z[i] ) )
	    {
	      if ( verbose ) std::cout << "Problem with file format at line " << 2 + atoms.names.size() << "."<< std::endl;
		  return false;
	    }
	  atoms.index[i] = i;
	  if (  !(same_type) || fill_lut ) addAtom2LUT( lut_list , { atoms.names[ i ] , atoms.x[ i ] , atoms.y[ i ] , atoms.z[ i ] , atoms.index[ i ] } );
	}
      i++;
    }
  //-------------------------------------------------------------------------------------


  return true; 
}
//---------------------------------------------------------------------------------------------------
bool readStepXYZ( std::ifstream & file , AtomList & atoms , AllTypeLUT & lut_list , const bool same_type , const bool verbose )
// Reads a step of XYZ file and build atom list and LUT for types
{
  //---------------------------
  // Clearing previous values
  //---------------------------------------------------------------
  if ( atoms.names.size() != 0 ) atoms.names.clear();
  if ( atoms.x.size() != 0 ) atoms.x.clear();
  if ( atoms.y.size() != 0 ) atoms.y.clear();
  if ( atoms.z.size() != 0 ) atoms.z.clear();
  if ( atoms.index.size() != 0 ) atoms.index.clear();
  if ( !(same_type) )
    {
      lut_list.types.clear();
      lut_list.type_index.clear();
      lut_list.type_name.clear();
    }
  else if ( lut_list.type_name.size() != 0 )
    {
      lut_list.types.clear();
      lut_list.type_index.clear();
      lut_list.type_name.clear();
    }
  //---------------------------------------------------------------
  
  //-------------
  // Filling LUT
  //---------------------------------------------------------------
  bool fill_lut;
  if ( lut_list.type_name.size() == 0 ) fill_lut=true;
  else fill_lut = false;
  //---------------------------------------------------------------

  //----------------------
  // Physical Parameters
  //---------------------------------------------------------------
  int nb_atoms = 0 ; // Numbers
  std::string line ; // Contains one line
  //---------------------------------------------------------------

  //------------------------------------------
  // Getting the number of atoms in the step
  //-----------------------------------------------------------------------------------
  if ( std::getline( file , line ) )
    {
      std::istringstream it_string(line);
      if ( ! ( it_string >> nb_atoms ) )
	{
	  if ( verbose ) std::cout << "Problem with file format at line 1." << std::endl;
	  return false;
	}
    }
  //-----------------------------------------------------------------------------------

  //-------------------------
  // Jumps next comment line
  //----------------------------------------------------------------------------------
  if( ! std::getline( file , line ) )
    {
      if ( verbose && line != "" )
	{
	  std::cout << "Problem with file format at line 2" << std::endl;
	}
      return false;
    }
  //----------------------------------------------------------------------------------

  //----------------------
  // Reads atoms by atoms
  //-------------------------------------------------------------------------------------
  while ( atoms.names.size() < nb_atoms )
    {
      if ( std::getline( file , line ) )
	{
	  std::istringstream it_string(line);
	  std::string name;
	  double x , y , z;
	  int index = atoms.names.size();
	  if ( ! ( it_string >> name >> x >> y >> z ) )
	    {
	      if ( verbose ) std::cout << "Problem with file format at line " << 2 + atoms.names.size() << "."<< std::endl;
	      return false;
	    }
	  atoms.names.push_back( name );
	  atoms.x.push_back( x ); atoms.y.push_back( y ); atoms.z.push_back( z );
	  atoms.index.push_back( index );
	  if ( !(same_type) || fill_lut )
	    {
	      addAtom2LUT( lut_list , { atoms.names[ index ] , atoms.x[ index ] , atoms.y[ index ] , atoms.z[ index ] , atoms.index[ index ] } );
	    }
	}
    }
  //-------------------------------------------------------------------------------------
  return true; 
}
//------------------------------------------------------------------------------------------------
bool readStepXYZ( std::ifstream & file , AtomList & atoms , std::vector<TypeLUT> & lut_list , const bool same_type , const bool verbose )
// Reads a step of XYZ file and build atom list and LUT for types
{
  //---------------------------
  // Clearing previous values
  //---------------------------------------------------------------
  if ( atoms.names.size() != 0 ) atoms.names.clear();
  if ( atoms.x.size() != 0 ) atoms.x.clear();
  if ( atoms.y.size() != 0 ) atoms.y.clear();
  if ( atoms.z.size() != 0 ) atoms.z.clear();
  if ( atoms.index.size() != 0 ) atoms.index.clear();
  if ( same_type && lut_list.size() != 0 ) lut_list.clear();
  bool fill_lut;
  if ( lut_list.size() == 0 ) fill_lut=true;
  else fill_lut = false;
  //---------------------------------------------------------------

  //----------------------
  // Physical Parameters
  //---------------------------------------------------------------
  int nb_atoms = 0 ; // Numbers
  Atom atom ;        // Stock a temporary atom
  std::string line ; // Contains one line
  //---------------------------------------------------------------

  //------------------------------------------
  // Getting the number of atoms in the step
  //-----------------------------------------------------------------------------------
  if ( std::getline( file , line ) )
    {
      std::istringstream it_string(line);
      if ( ! ( it_string >> nb_atoms ) )
	{
	  if ( verbose ) std::cout << "Problem with file format at line 1." << std::endl;
	  return false;
	}
    }
  //-----------------------------------------------------------------------------------

  //-------------------------
  // Jumps next comment line
  //----------------------------------------------------------------------------------
  if( ! std::getline( file , line ) )
    {
      if ( verbose && line != "" )
	{
	  std::cout << "Problem with file format at line 2" << std::endl;
	}
      return false;
    }
  //----------------------------------------------------------------------------------

  //----------------------
  // Reads atoms by atoms
  //-------------------------------------------------------------------------------------
  while ( atoms.names.size() < nb_atoms )
    {
      if ( std::getline( file , line ) )
	{
	  std::istringstream it_string(line);
	  std::string name;
	  double x,y,z;
	  int index = atoms.names.size();
	  if ( ! ( it_string >> name >> x >> y >> z ) )
	    {
	      if ( verbose ) std::cout << "Problem with file format at line " << 2 + atoms.names.size() << "."<< std::endl;
	      return false;
	    }
	  atoms.names.push_back( name );
	  atoms.x.push_back( x );
	  atoms.y.push_back( y );
	  atoms.z.push_back( z );
	  atoms.index.push_back( index );
	  if ( !(same_type) || fill_lut )
	    {
	      std::cout << "check" << std::endl;
	      addAtom2LUT( lut_list , name , index  );
	    }
	}
    }
  //-------------------------------------------------------------------------------------
  return true; 
}
//-------------------------------------------------------------------------------------------------
bool readStepXYZ( std::ifstream & file , std::vector<Atom> & atoms , std::vector<TypeLUT> & lut_list , const bool same_type , const bool verbose )
// Reads a step of XYZ file and build atom list and LUT for types
{

  // Clearing previous values
  //---------------------------------------------------------------
  if ( atoms.size() != 0 ) atoms.clear();
  if ( same_type && lut_list.size() != 0 ) lut_list.clear();
  //---------------------------------------------------------------

  // Physical Parameters
  //---------------------------------------------------------------
  int nb_atoms = 0 ; // Numbers
  Atom atom ;        // Stock a temporary atom
  std::string line ; // Contains one line
  //---------------------------------------------------------------

  // Getting the number of atoms in the step
  //-----------------------------------------------------------------------------------
  if ( std::getline( file , line ) )
    {
      std::istringstream it_string(line);
      if ( ! ( it_string >> nb_atoms ) )
	{
	  if ( verbose ) std::cout << "Problem with file format at line 1." << std::endl;
	  return false;
	}
    }
  //-----------------------------------------------------------------------------------

  // Jumps next comment line
  //----------------------------------------------------------------------------------
  if( ! std::getline( file , line ) )
    {
      if ( verbose && line != "" )
	{
	  std::cout << "Problem with file format at line 2" << std::endl;
	}
      return false;
    }
  //----------------------------------------------------------------------------------

  // Reads atoms by atoms
  //---------------------------------------------------------------
  while ( atoms.size() < nb_atoms )
    {
      if ( std::getline( file , line ) )
	{
	  std::istringstream it_string(line);
	  if ( ! ( it_string >> atom.name >> atom.x >> atom.y >> atom.z ) )
	    {
	      if ( verbose ) std::cout << "Problem with file format at line " << 2 + atoms.size() << "."<< std::endl;
	      return false;
	    }
	  atom.index = atoms.size();
	  atoms.push_back( atom );
	  if ( ! same_type ) addAtom2LUT( lut_list , atom );
	}
    }
  //---------------------------------------------------------------

  return true; 
}
//-----------------------------------------------------------------------------------
std::vector<Atom> readstepXYZ( std::ifstream & file )
{
  // Stream Handling
  std::istream_iterator<std::string> read(file);
  std::istream_iterator<std::string> end;
  
  // Atoms
  int nb_atoms=-1;
  Atom atom;
  std::vector<Atom> atom_list;
  int atom_count=0;
  
  // Reading step
  while( read != end && atom_list.size() < nb_atoms )
    {
      if ( nb_atoms > 0 )
	{
	  // Reads one atom
	  atom.name = std::string(*read); ++read;
	  atom.x    = atof(std::string(*read).c_str()); ++read;
	  atom.y    = atof(std::string(*read).c_str()); ++read;
	  atom.z    = atof(std::string(*read).c_str());
	  atom_list.push_back(atom);
	  if ( atom_list.size() != nb_atoms )
	    {
	      ++read;
	    }
	}
      else if ( nb_atoms == 0 )
	{
	  return atom_list;
	}
      else
	{
	  nb_atoms = atoi(std::string(*read).c_str());
	  ++read; // "STEP"
	  ++read; // STEP NUMBER
	  ++read;
	}
    }

  // Return atom list
  return atom_list;
}
//==================================================================================

//=======
// WRITE
//==================================================================================
void writeXYZ( std::ofstream & file , const std::vector<Atom> atom_list )
{
  file << atom_list.size() << std::endl;
  file << "STEP LOUTRE" << std::endl;
  for ( int i=0 ; i < atom_list.size() ; i++ )
    {
      file << atom_list[i].name << " " << atom_list[i].x <<  " " << atom_list[i].y <<  " " << atom_list[i].z << std::endl;
    }
  return;
}
//-----------------------------------------------------------------------------------
void writeXYZ( std::ofstream & file , const std::vector<Atom> atom_list , const int step )
{
  file << atom_list.size() << std::endl;
  file << "STEP " << step << std::endl;
  for ( int i=0 ; i < atom_list.size() ; i++ )
    {
      file << atom_list[i].name << " " << atom_list[i].x <<  " " << atom_list[i].y <<  " " << atom_list[i].z << std::endl;
    }
  return;
}
//-----------------------------------------------------------------------------------
void writeXYZ( std::ofstream & file , const AtomList atom_list )
{
  file << atom_list.x.size() << std::endl;
  file << "STEP X" << std::endl;
  for ( int i=0 ; i < atom_list.x.size() ; i++ )
    {
      file << atom_list.names[i] << " ";
      file << atom_list.x[i] << " ";
      file << atom_list.y[i] << " ";
      file << atom_list.z[i] << std::endl;
    }
  return;
}
//-----------------------------------------------------------------------------------
void writeXYZ( std::ofstream & file , const AtomList atom_list , const int step )
{
  file << atom_list.x.size() << std::endl;
  file << "STEP " << step << std::endl;
  for ( int i=0 ; i < atom_list.x.size() ; i++ )
    {
      file << atom_list.names[i] << " ";
      file << atom_list.x[i] << " ";
      file << atom_list.y[i] << " ";
      file << atom_list.z[i] << std::endl;
    }
  return;
}
//==================================================================================
