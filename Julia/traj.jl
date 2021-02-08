module traj

using utils
using periodicTable

export Traj

mutable struct Traj

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

    #---------------------------------------------------------------------------
    function Traj( nb_atoms::T1=0, nb_step::T2=0 ) where { T1 <: Int, T2 <: Int }
        dim=3
        if nb_atoms < 0
            print("Incorrect number of atoms. We have nb_atoms=", nb_atoms, " while it should be >= 0.\n" )
            return
        end
        if nb_step < 0
            print("Incorrect number of step. We have nb_step=", nb_step, " while it should be >= 0.\n" )
            return
        end
        new( Array{AbstractString,1}( undef, nb_atoms ),              # Atom_names
             Array{Int,1}(            undef, nb_atoms ),              # Atom_index
             Array{Real,3}(           undef, nb_step, nb_atoms, 3 ),  # Atom_position
             Array{Real,2}(           undef, nb_step, 3 ),            # Cell_lengths
             Array{Real,2}(           undef, nb_step, 3 ),            # Cell_angles
             Array{Real,3}(           undef, nb_step, 3, 3 ),         # Cell_matrix
            )
    end
    function Traj(  atom_name::T1, atom_index::T2, atom_positions::Vector{T3},
                    cell_lengths::Vector{T4}, cell_angles::Vector{T5}, cell_matrix::Array{T6,2} ) where { T1 <: AbstractString, T2 <: Int, T3 <: Real, T4 <: Real, T5 <: Real, T6 <: Real }
        # Description:
        # Creates a Traj structure containing a single atom and a single step
        # Input:
        # - atom_name (String): Name of the atom
        # - atom_index (Int): Index of the atom
        # - atom_positions ( Vector, length 3, Real): position of the atom in cell (default unit: Angstroms)
        # - cell_lengths (Vector, length 3, Real): a,b,c (lengths) of the sim cell (default unit: Angstroms)
        # - cell_angles (Vector, length 3, Real): alpha,beta,gamma (angles) of the sim cell (default unit: degrees)
        # - cell_matrix (Matrix, 3x3, Real): cell matrix of the sim cell (default unit: Angstroms)
        # Output:
        # - Trajectory with a single step and a single position
        positions = zeros(1,1,3)
        positions[1,1,:] = atom_positions
        cell_lengths_ = zeros(1,3)
        cell_lengths_ = cell_lengths
        cell_angles_  = zeros(1,3)
        cell_angles_  = cell_angles
        cell_matrix_  = zeros(1,3,3)
        cell_matrix_  = cell_matrix
        new( [atom_name], [atom_index], positions, cell_lengths_, cell_angles_, cell_matrix_ );
    end
    function Traj(  atom_name::T1, atom_index::T2, atom_positions::Vector{T3} ) where { T1 <: AbstractString, T2 <: Int, T3 <: Real }
        positions = zeros(1,1,3)
        positions[1,1,:] = atom_positions
        new( [atom_name], [atom_index], positions, zeros(1, 3), zeros(1,3), zeros(1,3,3) );
    end
    function Traj(  atom_name::T1, atom_positions::Vector{T3} ) where { T1 <: AbstractString, T2 <: Int, T3 <: Real }
        new( atom_name, zeros(1), atom_positions, zeros(1, 3), zeros(1,3), zeros(1,3,3) );
    end
    function Traj(  atom_names::Array{T1,1},
                    atom_index::Array{T2,1},
                    atom_positions::Array{T3,3},
                    cell_lengths::Array{T4,2},
                    cell_angles::Array{T5,2},
                    cell_matrix::Array{T6,3} ) where { T1 <: AbstractString, T2 <: Int, T3 <: Real, T4 <: Real, T5 <: Real, T6 <: Real }
        nb_atoms = size(atom_names)[1]
        nb_step = size(atom_positions)[1]
        if nb_atoms != size(atom_index)[1] || nb_atoms != size(atom_positions)[2]
            print( "Mismatch in number of atoms between data: \n" )
            print( "atom_names: ", nb_atoms,"\n" )
            print( "atom_index: ", size(atom_index)[1], "\n" )
            print( "atom_positions: ",size(atom_positions)[3] )
            print( "Aborting." )
            return
        end
        if nb_step != size(cell_lengths)[1] || nb_step != size(cell_angles)[1] || nb_step != size(cell_matrix)[1]
            print( "Mismatch in number of steps between data: \n" )
            print( "atom_positions: ", nb_step, "\n" )
            print( "cell_lengths: ", size(cell_lengths)[1], "\n" )
            print( "cell_angles: ", size(cell_angles)[1], "\n" )
            print( "cell_matrix: ", size(cell_matrix)[1], "\n" )
            print( "Aborting." )
            return
        end
        new( atom_names, atom_index, atom_positions, cell_lengths, cell_angles, cell_matrix );
    end
    #---------------------------------------------------------------------------
end

mutable struct TrajNVT

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

    function TrajNVT( atoms_names::Array{AbstractString,1}, atoms_index::Array{} )
    end

end

# Cell computation
#------------------------------------------------------------------------------
function cellmatrix2params( traj::T1, step::T2 )  where { T1 <: Traj, }
    for col=1:3
        for line=1:3
            traj.cell_length[col] += traj.cell_matrix[line,col]*traj.cell_matrix[line,col]
        end
        length[col] = sqrt( length[col] )
    end
    tau=180/pi
    angles = zeros( Real, 3 )
    angles[1] = acos(sum( cell_matrix[:,2].*cell_matrix[:,3] )/(length[2]*length[3]))*tau
    angles[2] = acos(sum( cell_matrix[:,1].*cell_matrix[:,3] )/(length[1]*length[3]))*tau
    angles[3] = acos(sum( cell_matrix[:,1].*cell_matrix[:,2] )/(length[1]*length[2]))*tau
    return Cell_param( length, angles )
end
#-------------------------------------------------------------------------------

end
