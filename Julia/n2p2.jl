module n2p2

# Description
# Functions used to deal with N2P2

export writeN2P2conf

function writeN2P2conf( positions::Array{T1,3}, forces::Array{T2,3}, cell_matrix::Array{T3,3}, atom_types::Vector{T4}, charge::Vector{T5}, charge_system::Vector{T6}, energy::Vector{T7}, file_out::T8 ) where { T1 <: Real, T2 <: Real, T3 <: Real, T4 <: AbstractString, T5 <: Real, T6 <: Real, T7<:Real, T8 <: AbstractString }
    # Arguments
    # - positions: positions of the atoms (nb_step,nb_atoms,3)
    # - forces: forces acting on atoms (nb_step)
    # - cell_matrix: cell tensor containing cell matrix with cell information
    # - atom_types: types of atoms
    # - charge: charges of atoms
    # - charge_system: total charge of the system
    # - energy: energy of the system
    # - file_out: file where to write out
    # Output
    # - Bool whether writting output worked

    # Get number of steps and atoms, and dimension
    nb_steps, nb_atoms, dim = size( positions )

    # Open output file
    handle_out = open(file_out,"w")

    # Loop over step
    for step=1:nb_steps
        # Write Begin signal
        write( handle_out, string( "begin\n" ) )

        # Write comment step
        write( handle_out, string( "comment STEP ", step, "\n" ) )

        # Loop over dimensions1
        for i=1:3
            # Write lattice signal
            write( handle_out, string( "lattice " ) )

            # Loop over dimensions2
            for j=1:3
                # Write cell matrix elements
                write( handle_out, string( cell_matrix[step,i,j], " " ) )
            end

            # Write end of line
            write( handle_out, string( "\n" ) )
        end

        # Loop over atoms
        for atom=1:nb_atoms
            # Write atom signal
            write( handle_out, string( "atom ")  )

            # Loop over dimensions
            for i=1:3
                # Write atomic positions
                write( handle_out, string( positions[step,atom,i], " " ) )
            end

            # Write atom type
            write( handle_out, string( atom_types[atom], " " ) )

            # Write atom charge
            write( handle_out, string( charge[atom], " " ) )

            # Write Unknown separator
            write( handle_out, string( "X " ) )

            # Loop over dimensions
            for i=1:3
                # Write forces on atom
                write( handle_out, string(forces[step,atom,i]," ") )
            end

            # Write end of line string
            write( handle_out, string( "\n" ) )
        end
        # Write energy
        write( handle_out, string( "energy ", energy[step], "\n") )

        # Write Charge of the system
        write( handle_out, string( "charge ", charge_system[step], "\n") )

        # Write end signal
        write( handle_out, string( "end\n" ) )
    end

    # Close file
    close(file)

    # Returns true if all went well to that point
    return true
end

end
