module metalwalls

# Load all necessary LibAtomicSim modules
using conversion
using LinearAlgebra
using atom_mod
using cell_mod
using filexyz
using pdb
using periodicTable
using utils

# Exporting functions
export writeData

function writeData( handle_out::T1,
    init_step::T2,
    nb_atoms::T3,
    cell::T4,
    atom_names::Vector{T5},
    atom_positions::Array{T6,2},
    nb_electrodes::T7=0 ) where { T1 <: IO, T2 <: Int, T3 <: Int, T4 <: cell_mod.Cell_param, T5 <: AbstractString, T6 <: Real, T7 <: Int }
    # Arguments
    # - handle_out: IO handler for output file
    # - init_step: number of the initial step
    # - nb_atoms: number of normal atoms
    # - cell : cell params with cell informations
    # - atom_names: vector (nb_atoms) with names of all (normal + electrodes) atoms
    # - atom_positions: array (3,nb_atoms) with positions of all (normal + electrodes) atoms
    # - nb_electrodes (opt): number of electrodes atoms
    # Output
    # - Bool: whether or not operation was successful

    # Write Header with system info
    write( handle_out, string("# header\n") )
    write( handle_out, string("step ", init_step, "\n") )
    write( handle_out, string("num_atoms ", nb_atoms, "\n") )
    write( handle_out, string("num_electrode_atoms ", nb_electrodes, "\n" ) )

    # Writes cell informations
    write( handle_out, string("#box\n") )
    for i=1:3
        write( handle_out, string(cell.lengths[i], " ") )
    end
    write( handle_out, string("\n") )

    # Write atomic positions
    write( handle_out, string("# coordinates\n") )
    for atom=1:size(atom_positions)[2]
        # Write atom name
        write( handle_out, string(atom_names[atom], " ") )

        # Loop over dimension
        for i=1:3
            # Write atomic coordinates, in Bohr
            write( handle_out, string(round( atom_positions[i,atom]*conversion.ang2Bohr, digits=3), " ") )
        end

        # Write end of line
        write( handle_out, string("\n") )
    end

    # Returns true if all went well
    return true
end
function writeData( file_path::T1,
    init_step::T2,
    nb_atoms::T3,
    cell::T4,
    atom_names::Vector{T5},
    atom_positions::Array{T6,2},
    nb_electrodes::T7=0 ) where { T1 <: AbstractString, T2 <: Int, T3 <: Int, T4 <: cell_mod.Cell_param, T5 <: AbstractString, T6 <: Real, T7 <: Int }
    # Arguments
    # - handle_out: path to the outpt file
    # - init_step: number of the initial step
    # - nb_atoms: number of normal atoms
    # - cell : cell params with cell informations
    # - atom_names: vector (nb_atoms) with names of all (normal + electrodes) atoms
    # - atom_positions: array (3,nb_atoms) with positions of all (normal + electrodes) atoms
    # - nb_electrodes (opt): number of electrodes atoms
    # Output
    # - Bool: whether or not operation was successful

    # Opening input file
    handle_out = open(file_path, "w")

    test = writeData( handle_out, init_step, nb_atoms, cell, atom_names, atom_positions, nb_electrodes=nb_electrodes )

    # Closing file
    close(handle_out)

    # Return status
    if test
        return true
    else
        return false
    end
end

function readRestart( handle_in::T1 ) where { T1 <: IO }
    # Argument
    # - handle_in: IO handler for input file
    # Output
    # - init_step: number of initial step before simulation
    # - nb_atoms: number of atoms
    # - nb_electrodes: number of electrode atoms
    # - atom_names: names of the atoms
    # - atom_positions: positions of the atoms
    # - cell : Cell_param with cell information

    # Skip first line
    readline( handle_in )

    # Read initial step
    init_step = parse(Int, split( readline( handle_in ) )[2] )

    # Read number of atoms
    nb_atoms = parse(Int, split( readline( handle_in ) )[2] )

    # Read number of electrode atoms
    nb_electrodes = parse(Int, split( readline( handle_in ) )[2] )

    # Skip box line
    readline( handle_in )

    # Get cell informations
    keys = split( readline( handle_in ) )

    # Initialize vectors for cell lengths
    lengths = zeros(Real, 3)

    # Loop over dimensions
    for i=1:3
        lengths[i] = parse(Float64, keys[i] )
    end

    # Converts to Cell_param
    cell = cell_mod.Cell_Param( lengths )

    # Initialize output vectors
    atom_names     = Vector{AbstractString}(undef, nb_atoms + nb_electrodes )
    atom_positions = zeros(Real, nb_atoms + nb_electrodes )

    # Loop over number of atoms
    for atom=1:(nb_atoms+nb_electrodes)
        # Read and parse atom line
        keys = split( readline( handle_in ) )

        # Put atom names
        atom_names[ atom ] = keys[1]

        # Loop over dimensions
        for i=1:3
            # Assign atomic positions
            atom_positions[i,atom] = keys[ i+1 ]
        end
    end

    # Return all information
    return init_step, nb_atoms, nb_electrodes, atom_names, atom_positions, cell
end
function readRestart( file_path::T1 ) where { T1 <: AbstractString }
    # Argument
    # - handle_in: IO handler for input file
    # Output
    # - init_step: number of initial step before simulation
    # - nb_atoms: number of atoms
    # - nb_electrodes: number of electrode atoms
    # - atom_names: names of the atoms
    # - atom_positions: positions of the atoms
    # - cell : Cell_param with cell information

    # Check file existence
    if ! isfile( handle_in )
        return false
    end

    # Opens file
    handle_in = open( file_path )

    # Reads file
    init_step, nb_atoms, nb_electrodes, atom_names, atom_positions, cell = readRestart( handle_in )

    # closes file
    close( handle_in )

    # REturns output
    return init_step, nb_atoms, nb_electrodes, atom_names, atom_positions, cell
end

end
