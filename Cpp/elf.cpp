//---------------
// External Libs
//-------------------
#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>  
#include <stdlib.h>    
#include <math.h>
//---------------------

//----------------
// Internal Files
//-------------------------
#include "utils.h"
#include "atom.h"
#include "cell.h"
#include "contact_matrix.h"
#include "xyz.h"
#include "histogram.h"
#include "pdb.h"
#include "sim.h"
#include "lut.h"
#include "molecules.h"
//-------------------------

//================
// MAIN PROGRAM
//=====================================================================
int main( void )
{
  //--------
  // Input
  //---------------------------------
  std::ifstream input("ELF.cube");
  //--------------------------------

  //--------
  // Output
  //-------------------------------------

  //-------------------------------------
  
  //---------------
  // Initializers
  //-------------------------------------
  AtomList  atom_list; // Atoms in cell
  //-------------------------------------

  //--------------------
  // Reading Cell File
  //-------------------------------------------------------------------
  Cell cell;
  if ( ! readParamCell( "cell.param" , cell ) )
    {
      return 1;
    }
  //-------------------------------------------------------------------

  std::string line;
  int nb_atoms=0;
  double ox, oy, oz;
  int nb_vox[3];
  double x,y,z;
  int atom_type;
  double dummy;
  double cell_x[3], cell_y[3], cell_z[3];
  AtomList atoms;
  
  //-------------------------------
  // Read useless two first lines
  //-------------------------------------------------------------------
  for ( int i=0; i<2; i++)
    {
  
      if ( ! std::getline(input,line) )
	{
	  std::cout << "Problem with the file ELF.cube" <<std::endl;
	  return 1;
	}
    }
  //-------------------------------------------------------------------

  //----------------------------------------------------------------------
  // Reading number of atoms and origin of the origin point of the density
  //------------------------------------------------------------------------
  if( std::getline(input,line) )
    {
      std::istringstream it_string(line);
      if ( ! ( it_string >> nb_atoms >> ox >> oy >> oz ) ) return 1;
    }
  else return 1;
  //-----------------------------------------------------------------------

  //------------------------------------------------------
  // Getting the index of the atoms to calculate the plots
  //--------------------------------------------------------
  int atom_1=-1;
  int atom_2=-1;
  do
    {
      std::cout << "Index of atom 1: " ;
      std::cin >> atom_1; 
      std::cout << "Index of atom 2: " ;
      std::cin >> atom_2;
    }
  while( atom_1 < 0 || atom_2 < 0 );
  //--------------------------------------------------------

  //-------------
  // Output file
  //----------------------------------------------------------------------------------
  // File name
  std::string output_name="output_"+(std::string)(atom_1)+"_"+(std::string)(atom_2)+".dat";
  // Opening file
  std::ofstream output(output_name.c_str());
  //----------------------------------------------------------------------------------

  
  //--------------------------------------------
  // Reading number of voxels and cell vectors 
  //-----------------------------------------------------------------------
  for( int i=0 ; i<3 ; i++ )
    {
      if( std::getline(input,line) )
	{
	  if ( ! ( it_string >> nb_vox[i] >> cell_x[i] >> cell_y[i] >> cell_z[i] ) ) return 1;
	}
      else return 1;
    }
  //-----------------------------------------------------------------------

  //----------------------------
  // Computing cell parameters
  //-----------------------------------------------------------------------
  Cell cell,cell2;
  cell.a = 0; cell.b=0; cell.c=0;
  for ( int i=0; i<3 ; i++ )
    {
      cell.a += cell_x[i]*cell_x[i];
      cell.b += cell_y[i]*cell_y[i];
      cell.c += cell_z[i]*cell_z[i];
    }
  cell.a = sqrt(cell.a)*nb_vox[0]*0.52917721092;
  cell.b = sqrt(cell.b)*nb_vox[1]*0.52917721092;
  cell.c = sqrt(cell.c)*nb_vox[2]*0.52917721092;
  cell2.a = sqrt(cell.a)*0.52917721092;
  cell2.b = sqrt(cell.b)*0.52917721092;
  cell2.c = sqrt(cell.c)*0.52917721092;
  //-----------------------------------------------------------------------
  
  //---------------------------
  // Reading atomic positions
  //----------------------------------------------------------------------------------------
  // Initiate the atom_list
  //-------------------------------------
  atoms.assign( nb_atoms, "");
  atoms.x.assign( nb_atoms, 0);
  atoms.y.assign( nb_atoms, 0);
  atoms.z.assign( nb_atoms, 0);
  atoms.index.assign( nb_atoms, 0);
  //-------------------------------------
  // Reading the atomic positions and all
  //---------------------------------------------------------------------------
  for( int i=0 ; i < nb_atoms ; i++ )
    {
      if ( std::getline(input,line) )
	{
	  if ( ! ( it_string >> atom_type >> dummy >> x >> y >> z ) ) return 1;
	}
      else return 1;
    }
  //-----------------------------------------------------------------------------


  // Reading Density/ELF
  //----------------------------------
  int element=0;
  while ( getline(input,line) )
    {
      if ()
    } 
  //----------------------------------
  
  return 0;
}
