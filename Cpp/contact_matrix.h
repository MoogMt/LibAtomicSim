#ifndef CONTACT_MATRIX_H
#define CONTACT_MATRIX_H

#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>  
#include <stdlib.h>    
#include <math.h>

#include "utils.h"
#include "atom.h"
#include "cell.h"
#include "cutoff.h"
#include "contact_matrix.h"

//================
// CONTACT MATRIX
//===================================================================================================================
// Restricted Matrix
//--------------------------------------------------------
struct Contact_Matrix
{
  std::vector<std::string> types;  // Types of the atoms
  std::vector<double> matrix;      // Contact matrix
};
//---------------------------------------------------------
// Full Contact Matrix
//--------------------------------------------------------
struct ContactMatrix
{
  int nb_atoms;               // Number of atoms
  AllTypeLUT lut_list;        // LUT list
  std::vector<double> matrix; // Contact_Matrix
};
//===================================================================================================================

//=======
// MAKE
//===================================================================================================================
// Restricted Matrix
//-------------------------------------------------------------------------------------------------
Contact_Matrix makeContactMatrix ( std::vector<Atom> & atom_list , const Cell box );
//------------------------------------------------------------------------------------------------
// Full Matrix
//-------------------------------------------------------------------------------------------------
// Connections - returns the contact matrix 
ContactMatrix makeContactMatrix ( AtomList & atom_list, const Cell cell , const CutOffMatrix cut_off , const AllTypeLUT lut_type );
// Distances - returns the contact_matrix
ContactMatrix makeContactMatrixDistance ( AtomList & atom_list, const Cell cell , const CutOffMatrix cut_off , const AllTypeLUT lut_type );
// Connections - modifies the cm
void makeContactMatrix ( ContactMatrix & cm , AtomList & atom_list, const Cell cell , const CutOffMatrix cut_off , const AllTypeLUT lut_list , const bool go_on = true );
// SoftConnection
ContactMatrix makeContactMatrixSoft ( AtomList & atom_list , const Cell cell , const CutOffMatrix cut_off , const AllTypeLUT lut_type , double r0 , int n, int m);
// Distances - modifies the cm
void makeContactMatrixDistance ( ContactMatrix & cm , AtomList & atom_list, const Cell cell , const CutOffMatrix cut_off , const AllTypeLUT lut_list , const bool go_on = true );
// Distances and Connection - modifies the contact matrices
void makeContactMatrix ( ContactMatrix & cm_connect , ContactMatrix & cm_distance , AtomList & atom_list, const Cell cell , const CutOffMatrix cut_off , const AllTypeLUT lut_list , const bool go_on = true );
//===================================================================================================================

//===================
// EXTRACTING MATRIX
//===================================================================================================================
ContactMatrix extractContactMatrix( const ContactMatrix old_cm , const std::string specie1 , const std::string specie2 );
ContactMatrix extractContactMatrix( const ContactMatrix old_cm , const std::vector<int> atom_list1 , const std::vector<int> atom_list2 );
//===================================================================================================================

//=============
// CONNECTION
//=====================================================================
bool connected( const ContactMatrix & cm , const int i , const int j);
//=====================================================================

//=========
// SORTING
//=============================================
/* sortContactMatrix
Aim: 
Input: A full contact matrix (ContactMatrix)
Output: Nothing
*/
void sortContactMatrix( ContactMatrix & cm );
//=============================================

//============
// DISTANCES
//=================================================================================================
// Restricted Contact Matrix
double getDistance( const Contact_Matrix & cm , const int atom_index1, const int atom_index2 );
//------------------------------------------------------------------------------------------------
// Full Contact Matrix
double getDistance( const ContactMatrix & cm , const int atom_index1 , const int atom_index2 );
//=================================================================================================

//========
// Angles
//===================================================================================================================
double getAngle( const Contact_Matrix & cm , const int atom_center_index , const int atom_2_index, const int atom_3_index );
//-------------------------------------------------------------------------------------------------------------------
double getAngle( const ContactMatrix & cm , const int atom_A , const int atom_B , const int atom_C );
//===================================================================================================================

//=========
// CONTACT
//===================================================================================================================
std::vector<double> getAtomContact( const Contact_Matrix & cm , const int atom_index );
std::vector<double> getAtomContact( const Contact_Matrix & cm , const int atom_index, const std::string specie );
void writeAtomContact( std::ofstream & file , const Contact_Matrix & cm , const std::vector<int> atom_index );
//-------------------------------------------------------------------------------------------------------------------
std::vector<double> getAtomContact( const ContactMatrix & cm , const int atom_index );
std::vector<double> getAtomContact( const ContactMatrix & cm , const int atom_index , const std::string specie );
//===================================================================================================================

//=============
// COORDINANCE
//===================================================================================================================
int getAtomNeighboursNb( const Contact_Matrix & cm , const int atom_index , const double cut_off_radius );
int getAtomNeighboursNb( const Contact_Matrix & cm , const int atom_index , const std::string specie , const double cut_off_radius );
std::vector<int> getAtomsNeighboursNb( const Contact_Matrix & cm , const std::vector<int> atom_index_list , const double cut_off_radius );
std::vector<int> getTypeNeighboursNb( const Contact_Matrix & cm , const std::string type , const double cut_off_radius ) ;
double getTypeCoordinance( const Contact_Matrix & cm , const std::string type , const double cut_off_radius );
//===================================================================================================================

//====================
// NEAREST NEIGHBORS
//==========================================================================================
// By index
double getNNearest( Contact_Matrix & cm , int n_nearest, int atom_index );
std::vector<double> getNNearest( Contact_Matrix & cm , std::vector<int> n_nearest, int atom_index );
std::vector<double> getNNearest( Contact_Matrix & cm , int nearest, std::vector<int> atom_indexes );
// By Specie
double getNNearest( Contact_Matrix & cm , int n_nearest, int atom_index, std::string specie );
std::vector<double> getNNearest( Contact_Matrix & cm , int n_nearest , std::string specie );
std::vector<double> getNNearest( Contact_Matrix & cm , int nearest, std::vector<int> atom_indexes , std::string specie );
std::vector<double> getNNearest( Contact_Matrix & cm , int n_nearest, std::vector<std::string> atom_types );
//------------------------------------------------------------------------------------------
// By Index
double getNNearest( ContactMatrix & cm , int n_nearest , int atom_index );
std::vector<double> getNNearest( ContactMatrix & cm , std::vector<int> n_nearest , int atom_index );
std::vector<double> getNNearest( Contact_Matrix & cm , int nearest , std::vector<int> atom_indexes );
std::vector<double> getNNearest( ContactMatrix & cm , int nearest, std::vector<int> atoms_center_index , std::vector<int> atoms_ext_index );
// By specie
double getNNearest( ContactMatrix & cm , int n_nearest, int atom_index, std::string specie );
std::vector<double> getNNearestVector( ContactMatrix & cm , int atom_index , int n_nearest , std::string specie );
std::vector<double> getNNearest( ContactMatrix & cm , int nearest, std::vector<int> atom_indexes , std::string specie );
std::vector<double> getNNearest( ContactMatrix & cm , int n_nearest, std::vector<std::string> atom_types );
//==========================================================================================

//=======
// IO
//=====================================================
// PRINT
//-----------------------------------------------------------------------------------------
void printContactMatrix( ContactMatrix & cm );
//-----------------------------------------------------------------------------------------
// WRITE
//-----------------------------------------------------------------------------------------
void writeAtomDistances( std::ofstream & file , std::vector<Atom> atom_list , std::vector<int> atom_index, Cell box );
void writeNearest( std::ofstream & file , Contact_Matrix contact_matrix , int n_nearest , std::vector<int> atom_indexes , int step );
void writeNearest( std::ofstream & file , Contact_Matrix contact_matrix , std::vector<int> nearest, int atom_index);
void writeNearest( std::ofstream & file , Contact_Matrix contact_matrix , std::vector<int> nearest, std::vector<int> atom_indexes , int step);
//==========================================================================================

#endif 
