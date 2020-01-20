#include "molecules.h"

//================//
// MAKE MOLECULES //
//=================================================================================
MoleculeBasic startMoleculeBasic( const std::string name_atom , const int index_atom )
{
  std::vector<std::string> names;
  std::vector<int> index;
  names.push_back( name_atom );
  index.push_back( index_atom );
  return { names, index };
}
//----------------------------------------------------------------------------
Molecule startMolecule( const std::string name_atom , const int index_atom )
{
  std::vector<std::string> names;
  std::vector<int> index;
  std::vector<Bond> bonds;
  names.push_back( name_atom );
  index.push_back( index_atom );
  return { names, index, bonds };
}
//----------------------------------------------------------------------------
std::vector<MoleculeBasic> makeMoleculesBasic( const ContactMatrix & cm )
{
  //-----------------
  // Initialization
  //------------------------------------------------------------------
  std::vector<MoleculeBasic> mol_list;
  // Initiate all connection vectors
  int used[ cm.nb_atoms ]; zeros( used , cm.nb_atoms ); // check if atom is already in a molecule
  int try1[ cm.nb_atoms ]; zeros( try1 , cm.nb_atoms ); // global connection vector
  int try2[ cm.nb_atoms ]; zeros( try2 , cm.nb_atoms ); // local atomic connection vector
  //------------------------------------------------------------------
  
  //-----------------------------------------------------------------------
  // Loop over all non already used atoms
  for ( int i=0 ; i < cm.nb_atoms ; i++ )
    {
      // if atom is already in a molecule, skip to the next one
      if ( used[i] == 1 ) continue;
      // Initialize connection matrix for atom i
      zeros( try1 , cm.nb_atoms ); try1[i] = 1;
      // Initiate molecule with atom i
      MoleculeBasic molecule = startMoleculeBasic( cm.lut_list.type_name[i] , i );
      // Looping to get all atoms in molecule
      do
	{
	  // Initiate connection vector for atom k
	  zeros( try2 , cm.nb_atoms ) ;
	  // Looping over all atoms in the connection vector
	  for ( int k=0 ; k < cm.nb_atoms ; k++ )
	    {
	      // If atom k isn't connected, skip it, else...
	      if ( try1[k] == 0 ) continue;
	      // If it is, mark it as used...
	      used[k] = 1;
	      // ... And look for its neighbours
	      for ( int h=0 ; h < cm.nb_atoms ; h++ )
		{
		  // If atom isn't used and is connected to local atom...
		  if ( used[h] == 0 && connected(cm,k,h) )
		    {
		      // Add it to the connection vector
		      try2[h] = 1;
		      // Add atom name to the molecule
		      molecule.names.push_back( cm.lut_list.type_name[h] );
		      // Add Atom index to the molecule
		      molecule.atom_index.push_back( h );
		    }
		}
	    }
	  // copy the new connection 
	  copy( try2 , try1 , cm.nb_atoms );
	  // If the connection doesn't give anything,
	  // then we have found all molecules.
	} while ( sum( try1, cm.nb_atoms ) != 0 );
      // Add molecule to list
      /*      if ( molecule.names.size() > 96 )
	{
	  std::cout << "Problem " << std::endl;
	  std::cout << "size1: " << molecule.names.size() << std::endl;
	  std::cout << "size2: " << molecule.atom_index.size() << std::endl;
	  std::cout << "index: " << molecule.atom_index[i] << std::endl;
	  }*/
      mol_list.push_back( molecule );
    }
  //-----------------------------------------------------------------------
  return mol_list;
}
//----------------------------------------------------------------------------
std::vector<Molecule> makeMolecules( const ContactMatrix & cm )
{
  //-----------------
  // Initialization
  //------------------------------------------------------------------
  std::vector<Molecule> mol_list;
  // Initiate all connection vectors
  int used[ cm.nb_atoms ]; zeros( used , cm.nb_atoms ); // check if atom is already in a molecule
  int try1[ cm.nb_atoms ]; zeros( try1 , cm.nb_atoms ); // global connection vector
  int try2[ cm.nb_atoms ]; zeros( try2 , cm.nb_atoms ); // local atomic connection vector
  //------------------------------------------------------------------
  
  //-----------------------------------------------------------------------
  // Loop over all non already used atoms
  for ( int i=0 ; i < cm.nb_atoms ; i++ )
    {
      // if atom is already in a molecule, skip to the next one
      if ( used[i] == 1 ) continue;
      // Initialize connection matrix for atom 
      zeros( try1 , cm.nb_atoms ); try1[i] = 1;
      // Initiate molecule with atom i
      Molecule molecule = startMolecule( cm.lut_list.type_name[i] , i );
      // Looping to get all atoms in molecule
      do
	{
	  // Initiate connection vector for atom k
	  zeros( try2 , cm.nb_atoms ) ;
	  // Looping over all atoms in the connection vector
	  for ( int k=0 ; k < cm.nb_atoms ; k++ )
	    {
	      // If atom k isn't connected, skip it, else...
	      if ( try1[k] == 0 ) continue;
	      // If it is, mark it as used...
	      used[k] = 1;
	      // ... And look for its neighbours
	      for ( int h=0 ; h < cm.nb_atoms ; h++ )
		{
		  // If atom isn't used and is connected to local atom...
		  if ( used[h] == 0 && connected(cm,k,h) )
		    {
		      // ... Add it to the connection vector
		      try2[h] = 1;
		      // ... And to the molecule
		      molecule.names.push_back( cm.lut_list.type_name[h] ); // Name
		      molecule.atom_index.push_back( h );                   // Index
		      Bond bond =  { k , h , 1.0 };
		      molecule.bonds.push_back( bond );                  // Bond
		    }
		}
	    }
	  // copy the new connection
	  copy( try2 , try1 , cm.nb_atoms );
	  // If the connection vector is null then we have found all atoms to the molecule
	} while ( sum( try1, cm.nb_atoms ) != 0 );
      // Add molecule to the molecule list
      /*std::cout << "index: " << molecule.atom_index.size() << std::endl; 
      std::cout << "size: " <<  molecule.names.size() << std::endl;
      for ( int i=0; i < molecule.atom_index.size()-1 ; i++ )
	{
	  for ( int j=i+1; j < molecule.atom_index.size() ; j++ )
	    {
	      if ( molecule.atom_index[i] == molecule.atom_index[j] )
		{
		  std::cout << "index_ : " <<  i << " " << molecule.atom_index[i] << std::endl;
		  std::cout << "index_ : " <<  j << " " << molecule.atom_index[j] << std::endl;
		}
	    }
	    }*/
      mol_list.push_back( molecule );
    }
  //-----------------------------------------------------------------------
  return mol_list;
}
//=================================================================================

//================
// PRINT MOLECULE
//=================================================================================
void printMolecule( const MoleculeBasic molecule )
{
  for ( int i=0 ; i < molecule.names.size() ; i++ )
    {
      std::cout << molecule.names[i] << " " << molecule.atom_index[i] << " ";
    }
  std::cout << std::endl;
  return;
}
//-----------------------------------------------------------------------------
void printMoleculeSize( const MoleculeBasic molecule, const bool toline )
{
  std::cout << molecule.names.size() << " ";
  if ( toline ) std::cout << std::endl;
  return;
}
//-----------------------------------------------------------------------------
void printMolecules( const std::vector<MoleculeBasic> molecules )
{
  std::cout << "------------------------" << std::endl;
  for ( int i=0 ; i < molecules.size() ; i++ )
    {
      printMolecule( molecules[i] );
    }
  std::cout << "------------------------" << std::endl;
  return ;
}
//-----------------------------------------------------------------------------
void printMoleculesSize( const std::vector<MoleculeBasic> molecules )
{
  for ( int i=0 ; i < molecules.size() ; i++ )
    {
      printMoleculeSize( molecules[i], false );
    }
  std::cout << std::endl;
  return;
}
//=================================================================================

//========
// Bonded
//=================================================================================
std::vector<int> getBonded( Molecule molecule , int atom_index )
{
  std::vector<int> bonded;
  for ( int i=0 ; i < molecule.bonds.size() ; i++ )
    {
      if ( molecule.bonds[i].atom1_index == atom_index ) bonded.push_back( molecule.bonds[i].atom2_index );
      else if ( molecule.bonds[i].atom2_index == atom_index ) bonded.push_back( molecule.bonds[i].atom1_index );  
    }
  return bonded;
}
//=================================================================================

//=========================
// GET ANGLES AROUND ATOMS
//=================================================================================
std::vector<double> getAngleAtom( ContactMatrix cm , Molecule molecule , int atom_index )
{
  std::vector<double> angles;
  std::vector<int> bonded_atom = getBonded( molecule , atom_index );
  if ( bonded_atom.size()  < 2 ) return angles;
  for ( int i=0 ; i < bonded_atom.size()-1 ; i++ )
    {
      for ( int j=i+1 ; j < bonded_atom.size() ; j++ )
	{
	  angles.push_back( getAngle( cm , atom_index , bonded_atom[i] , bonded_atom[j] ) );
	}
    }
  return angles;
}
//=================================================================================
