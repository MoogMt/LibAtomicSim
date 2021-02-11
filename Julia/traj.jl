module traj

using utils
using periodicTable
using atom_mod
using cell_mod

export Traj

# Description:
# Set of structures and associated functions to deal with trajectory, to be used
# for all input/ouput functions for Molecular dynamics softwares (GROMACS, LAMMPS,
# CPMD, CP2K, whatever) to manipulate data as structures of array instead of
# array of structures. 

# TODO: a whole lot, including
# - more structure for specific trajectory
# - more

# Structures
#----------------------------------------------------------
mutable struct Traj
    # Descroption
    # Structure that contains simultaneously all information in  array form
    # which allows for easier manipulation of trajectory by various functions

    # Variables
    #------------------------------------------------------
    # Atom information
    atoms_names        :: Array{AbstractString,1}
    atoms_index        :: Array{Int,1}
    atoms_positions    :: Array{Real,3}
    # Cell Information
    cell_lengths       :: Array{Real,2}
    cell_angles        :: Array{Real,2}
    cell_matrix        :: Array{Real,3}
    #------------------------------------------------------

    #--------------------------------------------------------------------------
    # Create a default traj, with a given amount of atoms and steps
    function Traj( nb_atoms::T1=0, nb_step::T2=0 ) where { T1 <: Int, T2 <: Int }
        # Arguments
        # - nb_atoms: number of atoms (int) ( optional, default=0 )
        # - nb_step: number of steps (int)  ( optional, default=0 )
        # Output
        # - Create a traj object of definite step and number of atoms or return false if nb_atoms or step is negative

        # Check that the given number of atoms is positive
        if nb_atoms < 0
            return false
        end

        # Check that the given number of steps is positive
        if nb_step < 0
            return false
        end

        # return the default Traj structure
        new( Array{AbstractString,1}( undef, nb_atoms ),              # Atom_names
             Array{Int,1}(            undef, nb_atoms ),              # Atom_index
             Array{Real,3}(           undef, nb_step, nb_atoms, 3 ),  # Atom_position
             Array{Real,2}(           undef, nb_step, 3 ),            # Cell_lengths
             Array{Real,2}(           undef, nb_step, 3 ),            # Cell_angles
             Array{Real,3}(           undef, nb_step, 3, 3 ),         # Cell_matrix
            )
    end
    # Creates a Traj with a single atom using its name, index and position (vector), with no definite cell, defining a default cell instead (lengths=1A, angles=90°)
    function Traj(  atom_name::T1, atom_index::T2, atom_positions::Vector{T3} ) where { T1 <: AbstractString, T2 <: Int, T3 <: Real }
        # Arguments
        # atom_name: name of the atom (string)
        # atom_index: index of the atom (int)
        # atom_positions: position of the atom (real vector of size 3)
        # Output
        # - Creates a trajectory object for the single atom

        # Converting the vector into a tensor form
        positions = zeros(1,1,3)
        positions[1,1,:] = atom_positions

        # Creating cell_matrix
        cell_matrix = zeros(1,3,3)
        cell_matrix[1,:,:] = Matrix{Real}(I,3,3)*1.0

        # Creating the traj object
        new( [atom_name], [atom_index], positions, ones(1, 3), ones(1,3)*90, cell_matrix );
    end
    # Creates a Traj witha  single atom using its name and positions (vector real), with no definite cell
    function Traj(  atom_name::T1, atom_positions::Vector{T3} ) where { T1 <: AbstractString, T2 <: Int, T3 <: Real }
        # Argument
        # - atom_name: name of the atom (string)
        # - atom_positions: position of the atom (vector of real, size 3)
        # Output
        # - A Traj object with a single object

        # Creates a new Traj with a single atom and a single step
        new( atom_name, zeros(1), atom_positions )
    end
    # Creates a Traj with several atom using their name, index and position (real matrix), with no definite cell, defining a default cell instead (lengths=1A, angles=90°)
    function Traj(  atom_names::Vector{T1}, atom_indexes::Vector{T2}, atom_positions::Array{T3,2} ) where { T1 <: AbstractString, T2 <: Int, T3 <: Real }
        # Arguments
        # atom_names: names of the atoms (vector of string)
        # atom_indexes: indexes of the atom (vector of int)
        # atom_positions: positions of the atoms (real array (nb_atoms,size 3)
        # Output
        # - Creates a trajectory object for the single atom

        # Get number of atoms
        nb_atoms = size(positions)[1]

        # Converting the vector into a tensor form
        positions = zeros( 1, nb_atoms, 3 )
        positions[1,1,:] = atom_positions

        # Creating cell_matrix
        cell_matrix = zeros(1,3,3)
        cell_matrix[1,:,:] = Matrix{Real}(I,3,3)*1.0

        # Creating the traj object
        new( atom_names, atom_indexes, positions, ones(1, 3), ones(1,3)*90, cell_matrix )
    end
    # Creates a Traj with an AtomList, with no definite cell, defining a default cell instead (lengths=1A, angles=90°)
    function Traj(  atoms::T1 ) where { T1 <: atom_mod.AtomList }
        # Arguments
        # atoms: AtomList that describes the atomic structure
        # Output
        # - Creates a trajectory object for the single atom

        # Creating the traj object
        Traj( atoms.names, atom.index, atom.positions )
    end
    # Creates a Traj structure containing a single atom and a single step in a definite cell
    function Traj(  atom_name::T1, atom_index::T2, atom_positions::Vector{T3}, cell_lengths::Vector{T4}, cell_angles::Vector{T5}, cell_matrix::Array{T6,2} ) where { T1 <: AbstractString, T2 <: Int, T3 <: Real, T4 <: Real, T5 <: Real, T6 <: Real }
        # Arguments:
        # - atom_name (String): Name of the atom
        # - atom_index (Int): Index of the atom
        # - atom_positions ( Vector, length 3, Real): position of the atom in cell (default unit: Angstroms)
        # - cell_lengths (Vector, length 3, Real): a,b,c (lengths) of the sim cell (default unit: Angstroms)
        # - cell_angles (Vector, length 3, Real): alpha,beta,gamma (angles) of the sim cell (default unit: degrees)
        # - cell_matrix (Matrix, 3x3, Real): cell matrix of the sim cell (default unit: Angstroms)
        # Output:
        # - Trajectory with a single step and a single position

        # Initialize positions to tensor form
        positions = zeros(1,1,3)
        positions[1,1,:] = atom_positions

        # Converts cell lengths to matrix form
        cell_lengths_ = zeros(1,3)
        cell_lengths_ = cell_lengths

        # Converts cell angles to matrix form
        cell_angles_  = zeros(1,3)
        cell_angles_  = cell_angles

        # Converts cell matrix to tensor form
        cell_matrix_  = zeros(1,3,3)
        cell_matrix_  = cell_matrix

        new( [atom_name], [atom_index], positions, cell_lengths_, cell_angles_, cell_matrix_ );
    end
    # Creates a Traj containing a single atom (AtomList), a single step and a single cell (Cell_param)
    function Traj(  atoms::T1, cell_params::T2 ) where { T1 <: atom_mod.AtomList, T2 <: cell_mod.Cell_param }
        # Arguments:
        # - atoms: AtomList containing a single atom
        # - cell_params: Cell_param containing the cell
        # Output:
        # - Trajectory with a single step and a single position

        # Returns the Trajectory
        new( atoms.names[1], atoms_index[1], atoms.positions[1,:], cell_params.lengths, cell_params.angles, cell_mod.params2Matrix(cell_params) );
    end
    # Creates a Traj from atom_names (vector), index(vector), positions(tensor),  cell lengths (matrix real), cell angles (matrix angles), matrix (real tensor)
    function Traj(  atom_names::Array{T1,1}, atom_index::Array{T2,1}, atom_positions::Array{T3,3}, cell_lengths::Array{T4,2}, cell_angles::Array{T5,2}, cell_matrix::Array{T6,3} ) where { T1 <: AbstractString, T2 <: Int, T3 <: Real, T4 <: Real, T5 <: Real, T6 <: Real }
        # Argument
        # - atom_names: vector of string with the name of atoms
        # - atom_index: vector of int with index of atoms
        # - atom_positions: real tensor (nb_step,nb_atoms,3) with the positions of the atoms
        # - cell_lengths: real matrix (nb_step,3) with the lengths of the cell
        # - cell_angles: real matrix (nb_step,3) with the angles of the cell
        # - cell_matrix: real tensor (nb_step,3,3) with the cell matrix
        # Output:
        # Traj object or false, if dimension mismatch

        # Get the number of atoms in the traj
        nb_atoms = size(atom_names)[1]

        # Get the number of step of the traj
        nb_step = size(atom_positions)[1]

        # Checking that the various size of array concerning the number of atoms match
        if nb_atoms != size(atom_index)[1] || nb_atoms != size(atom_positions)[2]
            # If fails, send information about the size of the various arrays
            print( "Mismatch in number of atoms between data: \n" )
            print( "atom_names: ", nb_atoms,"\n" )
            print( "atom_index: ", size(atom_index)[1], "\n" )
            print( "atom_positions: ",size(atom_positions)[3] )
            print( "Aborting." )
            # And return false
            return false
        end

        # Checking that the various size of array concerning number of step match
        if nb_step != size(cell_lengths)[1] || nb_step != size(cell_angles)[1] || nb_step != size(cell_matrix)[1]
            # If fails, send information about the size of the various arrays
            print( "Mismatch in number of steps between data: \n" )
            print( "atom_positions: ", nb_step, "\n" )
            print( "cell_lengths: ", size(cell_lengths)[1], "\n" )
            print( "cell_angles: ", size(cell_angles)[1], "\n" )
            print( "cell_matrix: ", size(cell_matrix)[1], "\n" )
            print( "Aborting." )
            # And return false
            return false
        end

        # If everything checks out, return the Traj
        new( atom_names, atom_index, atom_positions, cell_lengths, cell_angles, cell_matrix );
    end
    #---------------------------------------------------------------------------
end
mutable struct TrajNVT
    # Description:
    # Useful structure when one wants to deal with NVT calculation, in which you don't need
    # to have cell information for every step.

    # Variables
    #------------------------------------------------------
    # Atom information
    atoms_names        :: Array{AbstractString,1}
    atoms_index        :: Array{Int,1}
    atoms_positions    :: Array{Real,3}
    # Cell Information
    cell_lengths       :: Array{Real,1}
    cell_angles        :: Array{Real,1}
    cell_matrix        :: Array{Real,2}
    #------------------------------------------------------
end
#----------------------------------------------------------

end
