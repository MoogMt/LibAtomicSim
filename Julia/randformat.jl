module randformat

# Import all import module from LibAtomSim
using atom_mod
using cell_mod

# Import necessary library from julia default repository
using LinearAlgebra

# V1 Format
#==============================================================================#
function writeV1( handle_out::T1, atom_names::Vector{T2}, atom_positions::Array{T3,2}, cell_matrix::Array{T4,2} )  where { T1 <: IO, T2 <: AbstractString, T3 <: Real, T4 <: Real }
    # Arguments
    # - handle_out: IO handler of the output file
    # - atom_names: names of the atoms
    # - atom_positions: positions of the atoms (cartesian), format: (3,nb_atoms)
    # - cell_matrix: cell matrix
    # Output
    # - Bool: whether writting was successful

    # Write comment line
    write( handle_out, string("Unit cell vectors: \n") )

    # Write cell matrix
    #----------------------------------------
    # Write vector A
    write( handle_out, string("va= ") )
    for i=1:3
        write( handle_out, string( round( cell_matrix[1,i], digits=3 ), " " ) )
    end
    write( handle_out, string("\n") )
    #----------------------
    write( handle_out, string("vb= ") )
    for i=1:3
        write( handle_out, string( round( cell_matrix[2,i], digits=3 ), " " ) )
    end
    write( handle_out, string("\n") )
    #----------------------
    write( handle_out, string("vc= ") )
    for i=1:3
        write( handle_out, string( round( cell_matrix[3,i], digits=3 ), " " ) )
    end
    write( handle_out, string("\n") )
    #----------------------------------------

    # Write number of atoms
    write( handle_out, string( size( atom_names )[1], "\n" ) )

    # Loop over atoms
    for atom=1:size(atom_names)[1]
        # Write atom name
        write( handle_out, string( atom_names[atom], " ") )

        # Compute reduced coordinates
        position_reduced = cell_mod.cartesian2Reduced( atom_positions[:,atom], cell_matrix )

        # Loop over dimensions
        for i=1:3
            write( handle_out, string( round( position_reduced[i], digits=3 ), " " ) )
        end

        # Write end of line
        Base.write( handle_out, string("\n") )
    end

    # Returns true if all went well
    return true
end
function writeV1( file_path::T1, atom_names::Vector{T2}, atom_positions::Array{T3,2}, cell_matrix::Array{T4,2} )  where { T1 <: AbstractString, T2 <: AbstractString, T3 <: Real, T4 <: Real }
    # Arguments
    # - file_path: path to the output file
    # - atom_names: names of the atoms in the cell
    # - atom_positions: positions of the atoms format: (3,nb_atoms)
    # - cell_matrix: cell matrix
    # Output
    # - Bool: whether writting was successful

    # Opens file
    handle_out = open( file_path, "w" )

    # Writing file, gets status update
    test = writeV1( handle_out, atom_names, atom_positions, cell_matrix )

    # Closing file
    close( handle_out)

    # Returns true if all went well
    return test
end
function writeV1( handle_out::T1, atoms::T2, cell_matrix::Array{T3,2} )  where { T1 <: IO, T2 <: atom_mod.AtomList, T3 <: Real }
    # Arguments
    # - handle_out: IO handler of the output file
    # - atoms: AtomList with all atomic informations
    # - cell_matrix: cell matrix
    # Output
    # - Bool: whether writting was successful

    # Write comment line
    write( handle_out, string("Unit cell vectors: \n") )

    # Write cell matrix
    #----------------------------------------
    # Write vector A
    write( handle_out, string("va= ") )
    for i=1:3
        write( handle_out, string( round( cell_matrix[1,i], digits=3 ), " " ) )
    end
    write( handle_out, string("\n") )
    #----------------------
    write( handle_out, string("vb= ") )
    for i=1:3
        write( handle_out, string( round( cell_matrix[2,i], digits=3 ), " " ) )
    end
    write( handle_out, string("\n") )
    #----------------------
    write( handle_out, string("vc= ") )
    for i=1:3
        write( handle_out, string( round( cell_matrix[3,i], digits=3 ), " " ) )
    end
    write( handle_out, string("\n") )    #----------------------------------------

    # Write number of atoms
    write( handle_out, string( size( atoms.names )[1], "\n" ) )

    # Loop over atoms
    for atom=1:size(atoms.names)[1]
        # Write atom name
        write( handle_out, string( atoms.names[atom], " ") )

        # Compute reduced coordinates
        position_reduced = cell_mod.cartesian2Reduced( atoms.positions[:,atom], cell_matrix )

        # Loop over dimensions
        for i=1:3
            write( handle_out, string( round( position_reduced[i], digits=3 ), " " ) )
        end

        # Write end of line
        Base.write( handle_out, string("\n") )
    end

    # Returns true if all went well
    return true
end
function writeV1( file_path::T1, atoms::T2, cell_matrix::Array{T3,2} )  where { T1 <: AbstractString, T2 <: atom_mod.AtomList, T3 <: Real }
    # Arguments
    # - file_path: path to the output file
    # - atoms: AtomList with atomic informations
    # - cell_matrix: cell matrix
    # Output
    # - Bool: whether writting was successful

    # Opens file
    handle_out = open( file_path, "w" )

    # Writing file, gets status update
    test = writeV1( handle_out, atoms, cell_matrix )

    # Closing file
    close( handle_out)

    # Returns true if all went well
    return test
end
#==============================================================================#

end
