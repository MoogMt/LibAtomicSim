#include "contact_matrix.h"
//==================
// CONTACT MATRIX
//=======================================================================================
Contact_Matrix makeContactMatrix( std::vector<Atom> & atom_list , const Cell box)
// Constructs the restricted contact matrix
{
  Contact_Matrix contact_matrix;
  for ( int i=0 ; i < atom_list.size() ; i++ )
    {
      contact_matrix.types.push_back(atom_list[i].name);
      for (int j=i+1; j < atom_list.size() ; j++ )
	{
	  contact_matrix.matrix.push_back( distanceAtoms( atom_list , i , j , box) );
	}
    }
  return contact_matrix;
}
//--------------------------------------------------------------------------------------
ContactMatrix makeContactMatrix ( AtomList & atom_list , const Cell cell , const CutOffMatrix cut_off , const AllTypeLUT lut_type )
// Constructs the full contact matrix
{
  // Determines the number of atoms
  int nb_atoms = atom_list.x.size();

  // Initialize contact matrix with 0s
  std::vector<double> matrix; matrix.assign(nb_atoms*nb_atoms,0.);

  // Loop over all pairs of atoms
  for ( int i=0 ; i < nb_atoms-1 ; i++ )
    {
      for ( int j=i+1 ; j < nb_atoms ; j++ )
	{
	  // Getting the cut_off for the atom_i vs atom_j interaction
	  double cutoff = getCutOff( cut_off , lut_type.type_index[i] , lut_type.type_index[j] );
	  // Comparing distance to cut_off
	  if ( distanceAtomsSq( atom_list , i , j , cell) < cutoff*cutoff )
	    {
	      matrix[i*nb_atoms+j] = 1;
	      matrix[j*nb_atoms+i] = 1; 
	    }
	}
    }

  // Sending results
  return { nb_atoms, lut_type, matrix };
}
//--------------------------------------------------------------------------------------
ContactMatrix makeContactMatrixDistance ( AtomList & atom_list , const Cell cell , const CutOffMatrix cut_off , const AllTypeLUT lut_type )
// Constructs the full contact matrix
{
  // Determines the number of atoms
  int nb_atoms = atom_list.x.size();

  // Initialize contact matrix with 0s
  std::vector<double> matrix; matrix.assign(nb_atoms*nb_atoms,0.);

  // Loop over all pairs of atoms
  for ( int i=0 ; i < nb_atoms-1 ; i++ )
    {
      for ( int j=i+1 ; j < nb_atoms ; j++ )
	{
	  double dist=distanceAtomsSq( atom_list , i , j , cell);
	  matrix[i*nb_atoms+j] = dist;
	  matrix[j*nb_atoms+i] = dist; 
	}
    }

  // Sending results
  return { nb_atoms, lut_type, matrix };
}
//--------------------------------------------------------------------------------------
ContactMatrix makeContactMatrixSoft ( AtomList & atom_list , const Cell cell , const CutOffMatrix cut_off , const AllTypeLUT lut_type , double r0 , int n, int m )
// Constructs soft contact matrix
{
  // Determines the number of atoms
  int nb_atoms = atom_list.x.size();

  // Initialize contact matrix with 0s
  std::vector<double> matrix; matrix.assign(nb_atoms*nb_atoms,0.);

  // Loop over all pairs of atoms
  for ( int i=0 ; i < nb_atoms-1 ; i++ )
    {
      for ( int j=i+1 ; j < nb_atoms ; j++ )
	{
	  double dist = sigmoidPlumed( sqrt(distanceAtomsSq( atom_list , i , j , cell)) ,r0,n,m);
	  matrix[i*nb_atoms+j] = dist;
	  matrix[j*nb_atoms+i] = dist; 
	}
    }

  // Sending results
  return { nb_atoms, lut_type, matrix };
}
//--------------------------------------------------------------------------------------
void makeContactMatrix ( ContactMatrix & cm , AtomList & atom_list, const Cell cell , const CutOffMatrix cut_off , const AllTypeLUT lut_list , const bool go_on )
// Constructs the full contact matrix
{
  // Determines the number of atoms
  cm.nb_atoms = atom_list.x.size();

  // Initialize contact matrix with 0s
  if ( !(go_on) || cm.matrix.size() == 0 )
    {
      cm.matrix.assign(cm.nb_atoms*cm.nb_atoms,0.);
    }

  // Loop over all pairs of atoms
  for ( int i=0 ; i < cm.nb_atoms-1 ; i++ )
    {
      for ( int j=i+1 ; j < cm.nb_atoms ; j++ )
	{
	  // Computing cut off for i and j
	  double cutoff = getCutOff( cut_off , lut_list.type_index[i] , lut_list.type_index[j] );
	  // Comparing distance to cut_off
	  if ( distanceAtomsSq( atom_list , i , j , cell) < cutoff*cutoff )
	    {
	      cm.matrix[i*cm.nb_atoms+j] = 1;
	      cm.matrix[j*cm.nb_atoms+i] = 1; 
	    }
	  else
	    {
	      cm.matrix[i*cm.nb_atoms+j] = 0;
	      cm.matrix[j*cm.nb_atoms+i] = 0; 
	    }
	}
    }

  // Updating LUT
  cm.lut_list = lut_list;
  
  // Sending results
  return;
}
//--------------------------------------------------------------------------------------
void makeContactMatrixDistance ( ContactMatrix & cm , AtomList & atom_list, const Cell cell , const CutOffMatrix cut_off , const AllTypeLUT lut_list , const bool go_on )
// Constructs the full contact matrix
{
  // Determines the number of atoms
  cm.nb_atoms = atom_list.x.size();

  // Initialize contact matrix with 0so
  if ( !(go_on) || cm.matrix.size() == 0 )
    {
      cm.matrix.assign( cm.nb_atoms*cm.nb_atoms , 0. );
    }


  // Loop over all pairs of atoms
  for ( int i=0 ; i < cm.nb_atoms-1 ; i++ )
    {
      for ( int j=i+1 ; j < cm.nb_atoms ; j++ )
	{
	  double dist= sqrt( distanceAtomsSq( atom_list , i , j , cell) );
	  cm.matrix[i*cm.nb_atoms+j] = dist;
	  cm.matrix[j*cm.nb_atoms+i] = dist;
	}
    }

  // Updating LUT for types
  cm.lut_list = lut_list;

  return;
}
//--------------------------------------------------------------------------------------
void makeContactMatrix( ContactMatrix & cm_connect , ContactMatrix & cm_distance , AtomList & atom_list, const Cell cell , const CutOffMatrix cut_off , const AllTypeLUT lut_list , const bool go_on )
{
  // Determines the number of atoms
  cm_connect.nb_atoms = atom_list.x.size();
  cm_distance.nb_atoms = atom_list.x.size();

  // Initialize contact matrix with 0so
  if ( !(go_on) || cm_connect.matrix.size() == 0 )
    {
      cm_connect.matrix.assign( cm_connect.nb_atoms*cm_connect.nb_atoms , 0. );
    }
  if ( !(go_on) || cm_distance.matrix.size() == 0 )
    {
      cm_distance.matrix.assign( cm_distance.nb_atoms*cm_distance.nb_atoms , 0. );
    }

  // Loop over all pairs of atoms
  for ( int i=0 ; i < cm_connect.nb_atoms-1 ; i++ )
    {
      double offset1 = i*cm_distance.nb_atoms;
      double offset2 = i*cm_connect.nb_atoms;
      for ( int j=i+1 ; j < cm_connect.nb_atoms ; j++ )
	{
	  double distsq = distanceAtomsSq( atom_list , i , j , cell);
	  double dist= sqrt( distsq );
	  cm_distance.matrix[ offset1 + j ] = dist;
	  cm_distance.matrix[ j*cm_distance.nb_atoms + i ] = dist;
	  double cutoff = getCutOff( cut_off , lut_list.type_index[i] , lut_list.type_index[j] );
	  if ( distsq < cutoff*cutoff )
	    {
	      cm_connect.matrix[ offset2 + j ] = 1;
	      cm_connect.matrix[ j*cm_connect.nb_atoms + i ] = 1; 
	    }
	  else
	    {
	      cm_connect.matrix[ offset2 + j ] = 0;
	      cm_connect.matrix[ j*cm_connect.nb_atoms + i ] = 0; 
	    }
	}
    }

  // Updating LUT for types
  cm_connect.lut_list = lut_list;
  cm_distance.lut_list = lut_list;

  return;
}  
//=======================================================================================

//===================
// EXTRACTING MATRIX
//=================================================================================================
ContactMatrix extractContactMatrix( const ContactMatrix old_cm , const std::string specie1 , const std::string specie2 )
{
  ContactMatrix new_cm;
  new_cm.nb_atoms = 0;
  for ( int i=0 ; i < old_cm.nb_atoms ; i++ )
    {
      if ( old_cm.lut_list.type_name[i] != specie1 ) continue;
      int offset = old_cm.nb_atoms*i;
      for ( int j = 0 ; j < old_cm.nb_atoms ; j++ )
	{
	  if ( old_cm.lut_list.type_name[j] != specie2 ) continue;
	  else
	    {
	      new_cm.matrix.push_back( old_cm.matrix[ offset + j ] );
	      if ( j == 0 ) new_cm.nb_atoms++;
	    }
	}
    }
  return new_cm;
}
//-------------------------------------------------------------------------------------------------
ContactMatrix extractContactMatrix( const ContactMatrix old_cm , const std::vector<int> atom_list1 , const std::vector<int> atom_list2 )
{
  ContactMatrix new_cm;
  new_cm.nb_atoms = 0;
  for ( int i=0 ; i < atom_list1.size() ; i++ )
    {
      int offset = atom_list1[i]*old_cm.nb_atoms;
      for ( int j=0 ; j < atom_list2.size() ; j++ )
	{
	  new_cm.matrix.push_back( old_cm.matrix[ offset +j ] );
	  if ( j == 0 ) new_cm.nb_atoms++;
	}
    }
  return new_cm;
}
//=================================================================================================

//============
// GET SPECIE
//==================================================
std::vector<int> getSpecieIndex( const ContactMatrix & cm, std::string specie )
{
  std::vector<int> index;
  for ( int i=0; i < cm.lut_list.type_name.size() ; i++ )
    {
      if( cm.lut_list.type_name[i] == specie )
	{
	  index = cm.lut_list.types[i].atom_index;
	}
    }
  if ( index.size() == 0 )
    {
      std::cout << "Chemical specie " << specie << " was not found in the box." << std::endl; 
    }
  return index;
}
//==================================================

//===========
// CONNECTED
//==================================================
bool connected( const ContactMatrix & cm , const int i , const int j )
{
  if ( cm.matrix[ i*cm.nb_atoms + j ] ) return true;
  else return false;
}
//=================================================

//=================
// SORTING MATRIX
//=================================================
void sortContactMatrix( ContactMatrix & cm )
{
  for ( int i=0 ; i < cm.nb_atoms ; i ++ )
    {
      double offset = cm.nb_atoms*i;
      for ( int j=0 ; j < cm.nb_atoms-1 ; j++ )
	{
	  for ( int k=0 ; k < cm.nb_atoms ; k++ )
	    {
	      int index1 = offset + j, index2 = offset + k;
	      if ( cm.matrix[ index1 ] > cm.matrix[ index2 ] )
		{
		  switchV( cm.matrix , index1 , index2 );
		}
	    }
	}
    }
  return ;
}
//=================================================

//==========
// DISTANCE
//=======================================================================================
double getDistance( const Contact_Matrix & cm , const int atom_index_1 , const int atom_index_2 )
{
  int nb_atoms = cm.types.size();
  int mini=min(atom_index_1,atom_index_2);
  int maxi=max(atom_index_1,atom_index_2);
  return cm.matrix[computeSep(mini,nb_atoms)-mini+maxi-1];
}
//--------------------------------------------------------------------------------------
double getDistance( const ContactMatrix & cm , const int atom_index1 , const int atom_index2 )
{
  return cm.matrix[atom_index1*cm.nb_atoms+atom_index2];
}
//=======================================================================================

//==================
// ANGLES FUNCTIONS
//=======================================================================================
// Calculates the angles between three atoms centered on atom_center_index using contact matrix and Al-Kashi theorem
double getAngle( const Contact_Matrix & cm , const int atom_center_index , const int atom_2_index, const int atom_3_index )
// Get Angle
{
  double a = getDistance( cm , atom_center_index, atom_2_index);
  double b = getDistance( cm , atom_center_index, atom_3_index);
  double c = getDistance( cm , atom_2_index, atom_3_index);
  return acos( (a*a+b*b-c*c)/(2*a*b) )*180.0/M_PI;
}
//--------------------------------------------------------------------------------------
double getAngle( const ContactMatrix & cm , const int atom_A , const int atom_B , const int atom_C )
// Get Angle(A)=(BAC)
{
  double a = getDistance( cm , atom_A , atom_B );
  double b = getDistance( cm , atom_A , atom_C );
  double c = getDistance( cm , atom_B , atom_C );
  return acos( (a*a+b*b-c*c)/(2*a*b) )*180.0/M_PI;
}
//=======================================================================================

//================
// ATOMIC CONTACT
//=======================================================================================
std::vector<double> getAtomContact( const Contact_Matrix & cm , const int atom_index )
{
  std::vector<double> contact_atom;
  int nb_atoms = cm.types.size();
  int sep=computeSep(atom_index,nb_atoms);
  for ( int i=0 ; i < atom_index ; i++ )
    {
      contact_atom.push_back( cm.matrix[atom_index+sumBtw(nb_atoms-1,nb_atoms-i-1)-1]);
    }
  for ( int i=0; i < nb_atoms-atom_index-1 ; i++ )
    {
      contact_atom.push_back( cm.matrix[sep+i]);
    }
  return contact_atom;
}
//----------------------------------------------------------------------------------------
std::vector<double> getAtomContact( const Contact_Matrix & cm , const int atom_index, const std::string specie ) 
{
  std::vector<double> contact_atom;
  int nb_atoms = cm.types.size();
  for(int i=0 ; i < nb_atoms ; i++ )
    {
      if ( i != atom_index && cm.types[i] == specie )
	{
	  contact_atom.push_back( getDistance( cm , atom_index , i ) );
	}
    }
  return contact_atom ;
}
//-----------------------------------------------------------------------------------------
std::vector<double> getAtomContact( const ContactMatrix & cm , const int atom_index )
{
  std::vector<double> contact;
  int offset = cm.nb_atoms*atom_index;
  for ( int i=0 ; i < cm.nb_atoms ; i++ )
    {
      double element = cm.matrix[ offset + i ];
      if ( element > 0 ) contact.push_back( element );
    }
  return contact;
}
//-----------------------------------------------------------------------------------------
std::vector<double> getAtomContact( const ContactMatrix & cm , const int atom_index , const std::string specie )
{
  std::vector<double> contact;
  double offset = atom_index*cm.nb_atoms;
  for ( int i=0 ; i < cm.nb_atoms ; i++ )
    {
      double element = cm.matrix[ offset+i ];
      if ( cm.lut_list.type_name[i] == specie && element > 0 )
	{
	  contact.push_back( element );
	}
    }
  return contact;
}
//=======================================================================================

//=============
// COORDINANCE
//=======================================================================================
int getAtomNeighboursNb( const Contact_Matrix & cm , const int atom_index, const double cut_off_radius )
{
  int neighbours = 0;
  std::vector<double> contact = getAtomContact( cm , atom_index );
  for ( int i=0 ; i < contact.size() ; i++ )
    {
      if ( contact[i] < cut_off_radius )
	{
	  neighbours++;
	}
    }
  return neighbours;
}
//------------------------------------------------------------------------------------------------
int getAtomNeighboursNb( const Contact_Matrix & cm , const int atom_index, const std::string specie, const double cut_off_radius )
{
  int neighbours = 0;
  std::vector<double> contact = getAtomContact( cm, atom_index );
  for ( int i=0 ; i < contact.size() ; i++ )
    {
      if ( contact[i] < cut_off_radius && cm.types[i] == specie )
	{
	  neighbours++;
	}
    }
  return neighbours;
}
//------------------------------------------------------------------------------------------------
std::vector<int> getAtomsNeighboursNb( const Contact_Matrix & cm , const std::vector<int> atom_index_list , const double cut_off_radius )
{
  std::vector<int> neighbours_nb_list;
  for ( int i=0 ; i < atom_index_list.size() ; i++ )
    {
      neighbours_nb_list.push_back( getAtomNeighboursNb( cm , atom_index_list[i] , cut_off_radius ) );
    }
  return neighbours_nb_list;
}
//------------------------------------------------------------------------------------------------
std::vector<int> getTypeNeighboursNb( const Contact_Matrix & cm , const std::string type, const double cut_off_radius ) 
{
  int nb_atoms = cm.types.size();
  std::vector<int> neighbours_nb_list;
  for ( int i=0 ; i <  nb_atoms ; i++ )
    {
      if ( cm.types[i] == type )
	{
	  neighbours_nb_list.push_back( getAtomNeighboursNb( cm , i , cut_off_radius) );
	}
    }
  return neighbours_nb_list;
}
//------------------------------------------------------------------------------------------------
double getTypeCoordinance( const Contact_Matrix & cm , const std::string type, const double cut_off_radius )
{
  return average( getTypeNeighboursNb( cm , type , cut_off_radius ) );
}
//=======================================================================================

//===============================
// NEAREST NEIGHBOURS FUNCTIONS
//=======================================================================================
// RESTRICTED CONTACT MATRIX
//------------------------------------------------------------------------------------------------
double getNNearest( Contact_Matrix & cm , int n_nearest, int atom_index )
// Gets the n_nearest neighbours distance from atom with index atom_index using a contact matrix
{
  return sortVector( getAtomContact( cm , atom_index ) , true )[n_nearest-1];
}
//------------------------------------------------------------------------------------------------
std::vector<double> getNNearest( Contact_Matrix & cm , std::vector<int> n_nearest, int atom_index )
{
  std::vector<double> n_nearest_list;
  for( int i=0; i < n_nearest.size() ; i++ )
    {
      n_nearest_list.push_back( getNNearest( cm , n_nearest[i] , atom_index ) );
    }
  return n_nearest_list;
}
//------------------------------------------------------------------------------------------------
std::vector<double> getNNearest( Contact_Matrix & cm , int nearest, std::vector<int> atom_indexes )
{
  std::vector<double> n_nearest;
  for ( int i=0 ; i < atom_indexes.size() ; i++ )
    {
      n_nearest.push_back( getNNearest( cm , nearest, atom_indexes[i] ) );
    }
  return n_nearest;
}
//------------------------------------------------------------------------------------------------
// By Specie
//------------------------------------------------------------------------------------------------
double getNNearest( Contact_Matrix & cm , int n_nearest, int atom_index, std::string specie )
{
  return sortVector( getAtomContact( cm , atom_index, specie ) , true )[n_nearest-1];
}
//------------------------------------------------------------------------------------------------
std::vector<double> getNNearest( Contact_Matrix & cm , int n_nearest , std::string specie )
{
  return sortVector( getAtomContact( cm , n_nearest , specie ) , true );
}
//------------------------------------------------------------------------------------------------
std::vector<double> getNNearest( Contact_Matrix & cm , int n_nearest, std::vector<std::string> atom_types )
{
  std::vector<double> n_nearest_list;
  for ( int i=0 ; i < atom_types.size() ; i++ )
    {
      std::vector<double> temp = getNNearest( cm , n_nearest , atom_types[i] );
      for( int j=0 ; j < temp.size() ; j++ )
	{
	  n_nearest_list.push_back(temp[i]);
	}
    }
  return n_nearest_list;
}
//------------------------------------------------------------------------------------------------
std::vector<double> getNNearest( Contact_Matrix & cm , int nearest, std::vector<int> atom_indexes , std::string specie )
{
  std::vector<double> n_nearest;
  for ( int i=0 ; i < atom_indexes.size() ; i++ )
    {
      n_nearest.push_back( getNNearest( cm , nearest, atom_indexes[i] , specie) );
    }
  return n_nearest;
}
//------------------------------------------------------------------------------------------------
// FULL CONTACT MATRIX
//------------------------------------------------------------------------------------------------
double getNNearest( ContactMatrix & cm , int n_nearest , int atom_index )
{
  return sortVector( getAtomContact( cm , atom_index ) , true )[n_nearest-1];
}
//------------------------------------------------------------------------------------------------
std::vector<double> getNNearest( ContactMatrix & cm , std::vector<int> n_nearest, int atom_index)
{
  std::vector<double> n_nearest_list;
  for( int i=0; i < n_nearest.size() ; i++ )
    {
      n_nearest_list.push_back( getNNearest( cm, n_nearest[i] , atom_index ) );
    }
  return n_nearest_list;
}
//------------------------------------------------------------------------------------------------
std::vector<double> getNNearest( ContactMatrix & cm , int nearest, std::vector<int> atom_indexes)
{
  std::vector<double> n_nearest;
  for ( int i=0 ; i < atom_indexes.size() ; i++ )
    {
      n_nearest.push_back( getNNearest( cm , nearest , atom_indexes[i] ) );
    }
  return n_nearest;
}
//-----------------------------------------------------------------------------------------------
std::vector<double> getNNearest( ContactMatrix & cm , int nearest, std::vector<int> atoms_center_index , std::vector<int> atoms_ext_index )
{
  std::vector<double> n_nearest;
  // Loop over atoms at the center
  for ( int i=0 ; i <  atoms_center_index.size() ; i++ )
    {
      // Contact Restricted to atom of interest
      std::vector<double> restricted_contact;
      // Offset to speed up calculation
      double offset = atoms_center_index[i]*cm.nb_atoms;
      // Loop over peripheral atoms of interest
      for ( int j=0 ; j < atoms_ext_index.size() ; j++ )
	{
	  // Distance between atoms
	  double element = cm.matrix[ offset + atoms_ext_index[j] ];
	  // Updating restricted contact
	  if ( element > 0) restricted_contact.push_back( element );
	}
      // Add nearest_th distance to the distance vector
      n_nearest.push_back( sortVector( restricted_contact , true )[ nearest-1 ] );
    }
  // Return the vector
  return n_nearest;
}
//------------------------------------------------------------------------------------------------
// By Specie
//------------------------------------------------------------------------------------------------
double getNNearest( ContactMatrix & cm , int atom_index , int nearest , std::string specie )
{
  return sortVector( getAtomContact( cm , atom_index, specie ) , true )[ nearest-1 ];
}
//------------------------------------------------------------------------------------------------
std::vector<double> getNNearestVector( ContactMatrix & cm , int atom_index , int n_nearest , std::string specie )
{
  
  return sortVector( getAtomContact( cm , atom_index, specie ) , true );
}
//------------------------------------------------------------------------------------------------
std::vector<double> getNNearest( ContactMatrix & cm , int nearest, std::vector<int> atom_indexes , std::string specie )
{
  std::vector<double> n_nearest;
  for ( int i=0 ; i < atom_indexes.size() ; i++ )
    {
      n_nearest.push_back( getNNearest( cm ,  atom_indexes[i] , nearest , specie) );
    }
  return n_nearest;
}
//------------------------------------------------------------------------------------------------
std::vector<double> getNNearest( ContactMatrix & cm , int atom_index, int n_nearest, std::vector<std::string> atom_types )
{
  std::vector<double> n_nearest_list;
  for ( int i=0 ; i < atom_types.size() ; i++ )
    {
      std::vector<double> temp = getNNearestVector( cm , atom_index , n_nearest , atom_types[i] );
      for( int j=0 ; j < temp.size() ; j++ )
	{
	  n_nearest_list.push_back(temp[i]);
	}
    }
  return n_nearest_list;
}
//===================================================================================================

//======
// IO
//==================================================================================================
// PRINT
//------------------------------------------------------------------------------------------------
void printContactMatrix( ContactMatrix & cm )
{
  for (int i=0; i < cm.nb_atoms*cm.nb_atoms; i++)
    {
      std::cout << cm.matrix[i] << " ";
      if ( (i+1) % cm.nb_atoms == 0 )
	{
	  std::cout << std::endl;
	}
    }
}
//------------------------------------------------------------------------------------------------
// WRITE
//------------------------------------------------------------------------------------------------
void writeAtomContact( std::ofstream & file , Contact_Matrix contact_matrix , std::vector<int> atom_indexes )
// Writes the atomic contact of a set of atoms in file
{

  for ( int i=0 ; i < atom_indexes.size() ; i++ )
    {
      std::vector<double> contact = getAtomContact(contact_matrix,atom_indexes[i]);
      
      for ( int j=0 ; j < contact.size() ; j++ )
	{
	  file << contact[j] << " " ;
	}
    }
  file << std::endl;
  return;
}
//-----------------------------------------------------------------------------------------
void writeNearest( std::ofstream & file , Contact_Matrix contact_matrix , int n_nearest , std::vector<int> atom_indexes , int step )
{
  file << step << " ";
  for ( int i=0 ; i < atom_indexes.size() ; i++ )
    {
      file << getNNearest( contact_matrix , n_nearest , atom_indexes[i]) << " ";
    }
  file << std::endl;
  return;
}
//------------------------------------------------------------------------------------------------
void writeNearest( std::ofstream & file , Contact_Matrix contact_matrix , std::vector<int> n_nearest, int atom_index)
// Writes the n_nearest atoms of atoms with index atom_index in file using contact matrix
{
  for( int i=0; i < n_nearest.size() ; i++ )
    {
      file << getNNearest( contact_matrix , n_nearest[i] , atom_index ) << " " ;
    }
  return;
}
//------------------------------------------------------------------------------------------------
void writeNearest( std::ofstream & file , Contact_Matrix contact_matrix , std::vector<int> nearest, std::vector<int> atom_indexes , int step)
// Gets and writes n earest neighbours distances 
{
  file << step << " ";
  for ( int i=0 ; i < atom_indexes.size() ; i++ )
    {
      writeNearest( file , contact_matrix , nearest, atom_indexes[i] );
    }
  file << std::endl;
  return;
}
//==================================================================================================
