module quantumespresso

# Description
# - Functions that deals with input/output specific to quantum espresso

# Export functions
export getEnergy
export computeBindingEnergy
export getMagnetization

# Get the energy from a simulation
#------------------------------------------------------------------------------#
function getEnergy( file_path::T1 ) where { T1 <: AbstractString }
    # Argument
    # - file_path: path to the output file of Quantum Espresso
    # Output
    # - energy: energies of the simulation

    # Initialize energy list (as empty)
    energy=[]

    # Opens file
    handle_in = open( file_path, "r" )

    # Loop over all lines
    while ! eof( handle_in )
        # Read each line
        line = readline( handle_in )

        # Check if the line contains the "!" flag specific to energy
        if contains( line, "!" )
            # Parse line to get the energy (the idea is to get the element between = and R (for Ry) )
            push!( energy, parse(Float64, AbstractString( split( split( line, "=" )[2], "R" )[1 ]) ) )
        end
    end

    # Close input file
    close( handle_in )

    # Returns energy in eV
    return energy .* conversion.ry2eV
end
#------------------------------------------------------------------------------#

# get relax binding energy
#------------------------------------------------------------------------------#
function computeBindingEnergy( file_pathj::T1, total_energy::T2, specie_number::Vector{T3}, isolated_energy::Vector{T4} ) where { T1 <: AbstractString, T2 <: Real, T3 <: Int, T4 <: Real }
    # Arguments
    # - file_path: path of the output file of QE
    # - specie_number: number of elements of each species
    # - isolated_energy: isolated energy of each species
    # Output
    # - Binding energy of the system

    # Get number of species
    nb_species = size( specie_number )[1]

    # Check that the number of species matches the dimension of the number of isolated energies
    if nb_species != size(isolated_energy)[1]
        # If they do not match, return false and sends a message
        error("Vector sizes do not match. Stopping now.\n")
        return false
    end

    # Initialize number of atoms to 0
    atoms_nb=0

    # Loop over species
    for specie=1:nb_species
        # Compute reduced energy
        total_energy = total_energy - specie_number[i]*isolated_energy[i]

        # Compute total number of atoms
        atoms_nb = atoms_nb + specie_number[i]
    end

    # Returns biding energy of the system
    return -energy/atoms_nb*conversion.ry2eV
end
#------------------------------------------------------------------------------#

# Magnetization
#------------------------------------------------------------------------------#
function getMagnetization( file_path::T1 ) where { T1 <: AbstractString }
    # Argument
    # - file_path: path of the input file
    # Output
    # - mag: magnetization of the system

    mag = zeros(Real, 0)

    # Open input file
    handle_in = open( file_path, "r" )

    # Loop over existing files
    while ! eof( handle_in )
        # Read line over each line of the file
        line = readline( handle_in )

        # Check that the line contains "total_magnetization"
        if contains( line, "total magnetization" )
            # GEt the string between "=" and "Bohr", and converts it into a number (the magnetization)
            push!( mag, parse(Float64, AbstractString( split( split( line, "=" )[2], "Bohr" )[1] ) ) )
        end
    end

    # Close input file
    close( handle_in )

    # Return magnetization
    return mag
end
#------------------------------------------------------------------------------#

end
