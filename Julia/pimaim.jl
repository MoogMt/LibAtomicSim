module pimaim

# Load all necessary LibAtomicSim modules
using conversion
using LinearAlgebra
using atom_mod
using cell_mod
using filexyz
using pdb
using periodicTable
using utils

# Export the function in the file
export getSpeciesAndNumber
export readPositions, readPositionsUpToCrash
export readPosCar
export readXV
export readRestart
export readCellParams, readCellBox
# NB: check that all funtions are here

# Description
# Set of functions that deals with input and output of the polarizable force
# field software called PIMAIM created by Mathieu Salanne's team in PHENIX
# in Sorbonne University, Paris, France (contact him if you want to know more)

# TODO: Test functions to check what works and what needs works

# Reads runtime.inpt input to get the species and number of species
#-------------------------------------------------------------------------------
function getSpeciesAndNumber( path_file::T1 ) where { T1 <: AbstractString }
    # Argument:
    # - path_file: the path to runtime.inpt
    # Output:
    # - species: vector of string, containing the species names or false if something goes wrong
    # - species_nb: vector of int, containing the number of element for each chemical specie, or false is something goes wrong

    # Check if file exist, otherwise return two false
    if ! isfile( path_file )
        return false, false
    end

    # Open input file
    handle_in = open( path_file )

    # Skips the first 4 lines
    utils.skipLines( handle_in, 4 )

    # Compute number of species, by reading the 5th line
    nb_species=parse( Int, split( readline( handle_in ) )[1] )

    # Initialize a vector for the species names
    species=Vector{AbstractString}(undef,nb_species)

    # Initialize a vector for the species number
    species_nb = zeros(Int, nb_species )

    # Read the line containing species names
    species_line = split( readline( handle_in ), "," )

    # Loop over all the element except the last one
    for i_spec = 1:nb_species-1
        # Get the species names (except for the last one)
        species[i_spec] = species_line[i_spec]
    end
    # Get the last remaining specie's name
    species[ nb_species ] = split(species_line[nb_species])[1]  # Avoid pesky comment

    # Read the line containing the species number
    species_nb_line = split( readline( handle_in ),"," )

    # Loop over the species
    for i_spec = 1:nb_species-1
        # Get the number of element for each chemical element
        species_nb[i_spec] = parse(Int, species_nb_line[i_spec] )
    end
    # Get number of element for final specie
    species_nb[ nb_species ] = parse(Int, split(species_nb_line[nb_species])[1] )  # Avoid pesky comment

    # Close the file
    close( handle_in )

    # Return the species, and their number
    return species, species_nb
end
#-------------------------------------------------------------------------------

# Reading positions
#==============================================================================#

# Functions dealing with POSCAR file
#-----------------------------------------------------------------------------
# Reads POSCAR file, returns AtomList and cell matrix corresponding to the last step of the simulation
function readPOSCAR( path_file::T1 ) where { T1 <: AbstractString }
    # Argument:
    # - path_file: path to the POSCAR file
    # Output
    # - atoms: AtomList contaning all atomic coordinates
    # - matrix: cell matrix, in matrix form (3x3, real)
    # Or false, false if there is a problem with the file

    # Get the number of lines of the file
    nb_lines = utils.getNbLines( path_file )
    # If the file is empty or the file does not exist, return two false
    if nb_lines == false || nb_lines == 0
        return false, false
    end

    # Open the file
    handle_in = open( path_file )

    # Skip the first two lines
    utils.skipLines( handle_in, 2)

    # Initialize the cell matrix
    matrix=zeros(Real,3,3)

    # Loop over the next three lines
    for i=1:3
        # Parse the line
        keyword = split( readline( handle_in) )
        # Loop over first 3 elements
        for j=1:3
            # Construct the matrix with the elements
            matrix[i,j] = parse(Float64,keyword[j])
        end
    end

    # Makes a copy of the matrix
    matrix2  = copy( matrix )

    # Create the reduced matrix of pimaim
    # - Loop over dimensions
    for i=1:3
        # Reduce all vectors by their norm
        matrix2[i,:] /= LinearAlgebra.norm( matrix2[i,:] )
    end

    # NB: Unsure about this
    # Check that the reduced matrix is ok
    if matrix2[2,1] == 0
        # If the matrix is not oriented properly, get the transpose
        matrix2 = transpose( matrix2 )
    end

    # Read the next line, parse it to get the species in the structure
    species = split( readline( handle_in ) )

    # Read the next line and by parsing it get the number of element per species
    nb_species_lines = split( readline(handle_in ) )

    # Initialize a vector of int for number of element per specie
    nb_species = zeros( Int, size(species)[1] )

    # Loop over each specie
    for i_spec=1:size(species)[1]
        # Convert the string into int for each element per specie
        nb_species[i_spec] = parse( Int, nb_species_lines[i_spec] )
    end

    # Construct the vector of string for the names of atoms based on the number of element and their order
    names_ = atom_mod.buildNames( species, nb_species )

    # Skipping line
    readline( handle_in )

    # Compute number of atoms in the structure
    nb_atoms=sum(nb_species)

    # Initialize AtomList with a given number of atoms
    atoms=AtomList(nb_atoms)

    # Initialize a temporary vector for reduced positions
    temp = zeros(Real,3)

    # Loop over atoms
    for atom=1:nb_atoms
        # Read the line corresponding to each atom position
        keys = split( readline( handle_in ) )

        # Loop over dimension
        for i=1:3
            # Parse into Float the string
            temp[i] = parse( Float64, keys[i] )
        end

        # Loop over dimension
        for i=1:3
            # Loop over dimension for matrix product to convert reduced into cartesian coordinates
            for j=1:3
                # Transform reduced coordinates into cartesian
                atoms.positions[ i, atom ] += matrix2[ j, i ]*temp[j]
            end
        end

        # Affect the name of the atom to the AtomList section
        atoms.names[atom] = names_[atom]

        # Affect the index of the atom
        atoms.index[atom] = atom

    end

    # Close the file
    close( handle_in )

    # Return AtomList with atomic information and cell matrix for the cell
    return atoms, matrix
end
#-----------------------------------------------------------------------------

# Functions dealing with poscart.out
#-----------------------------------------------------------------------------
# Reads poscart.out and returns a traj in term of vector of AtomList, position may still be in pseudo-reduced form
function readPoscarOut( path_file::T1, species::Vector{T2}, nb_species::Vector{T3} ) where { T1 <: AbstractString, T2 <: AbstractString, T3 <: Int }
    # Argument
    # - path_file: path to the poscart.out file (string)
    # - species: names of the species in the cell (vector of string)
    # - nb_species: number of element for each species (vector of int)
    # Output
    # - traj: vector of AtomList of the trajectory or false if file does not exists
    # OR returns false if poscart.out does not exists

    # Get number of lines of the file
    nb_lines = utils.getNbLines( path_file )
    # If the file does not exists, or is empty, return false
    if nb_lines == false || nb_lines == 0
        return false
    end

    # Compute the number of atoms in the traj
    nb_atoms = sum(nb_species)

    # Construct the names of atoms using the names of species and number of element per specie
    names_ = atom_mod.buildNames( species, nb_species )

    # Check that the format is ok (basic check for sanity of file)
    if nb_lines % nb_atoms != 0
        # If the format is problematic, returns false
        print( "Problem within file: ", path_file, " :\n" )
        print( "Number of lines is not a multiple of the number of atoms given.\n" )
        return false
    end

    # Compute the number of step in the trajectory
    nb_step = Int(nb_lines/nb_atoms)

    # Initialize the vector of AtomList for the positions
    traj = Vector{AtomList}( undef, nb_step )

    # Open the file
    handle_in = open( path_file )

    # Loop over step
    for step=1:nb_step

        # Initialize the AtomList of the step
        traj[step] = atom_mod.AtomList(nb_atoms)

        # Loop over atoms
        for atom=1:nb_atoms

            # Read and parse each line per atom
            keys = split( readline( handle_in ) )

            # Loop over index
            for i=1:3
                # Parse the positions, converts strings to float and converts from bohr to angstroms
                traj[ step ].positions[ i, atom ] = parse( Float64, keys[i] )*conversion.bohr2Ang
            end

            # Put the name of the atoms in
            traj[ step ].names[ atom ] = names_[atom]

            # Put the index of the atoms in
            traj[ step ].index[ atom ] = atom
        end

    end

    # Close the file
    close( handle_in )

    # Return the trajectory in the form of a vector of AtomList, position may need to be transformed
    return traj
end
# Reads poscart.out and returns a traj in term of vector of AtomList, use cell_matrices to transform the positions
# NB: unsure of the purpose of this function
function readPoscarOut( path_file::T1, species::Vector{T2}, nb_species::Vector{T3}, cell_matrices::Array{T4,3} ) where { T1 <: AbstractString, T2 <: AbstractString, T3 <: Int, T4 <: Real }
    # Argument
    # - path_file: path to the poscart.out file (string)
    # - species: names of the species in the structure (vector of string)
    # - nb_species: number of atoms per elements (vector of int)
    # - cell_matrices: tensor contaning the information of the cell (real matrix, nb_step,3,3)
    # Output
    # - traj: vector of AtomList, containing the trajectory of the system

    # Get the number of lines of the file
    nb_lines = utils.getNbLines( path_file )
    # If the file is empty or does not exists, return false
    if nb_lines == false
        return false
    end

    # Computes the number of atoms in the traj
    nb_atoms = sum(nb_species)

    # Construct the names
    names_ = atom_mod.buildNames( species, nb_species )

    # Check that the file format is ok
    if nb_lines % nb_atoms != 0
        # If the format is problematic, sends a message and returns false
        print("Problem within file: ",path_file," :\n")
        print("Number of lines is not a multiple of the number of atoms given.\n")
        return false
    end

    # Computes number of step of trajectory
    nb_step = Int( nb_lines/nb_atoms )

    # Initialize trajectory
    traj = Vector{AtomList}( undef, nb_step )

    # Open poscart.out file
    handle_in = open( path_file )

    # Loop over step
    for step=1:nb_step

        # Initialize local AtomList
        traj[step] = atom_mod.AtomList(nb_atoms)

        # Is
        matrix = copy( cell_matrices[step,:,:] )

        # Invert the given matrix
        inv_matrix = inv( matrix )

        # Get the proper matrix to convert the positions
        matrix2 = cell_mod.params2Matrix( cell_mod.cellMatrix2Params( matrix ) )

        # Loop over atoms
        for atom=1:nb_atoms

            # Reads and parse the atomic positions lines
            keys = split( readline( handle_in ) )

            # Loop over dimension
            for i=1:3
                # Read, parse positions from string to real, and converts from bohr to angstrom
                traj[step].positions[ i, atom ] = parse( Float64, keys[i] )*conversion.bohr2Ang
            end

            # conversion of positions
            positions_temp = matrix2 * inv_matrix * traj[step].positions[ :, atom ]

            # Copies atomic positions to AtomList
            traj[ step ].positions[ :, atom ] = positions_temp

            # Copies names of atom to AtomList
            traj[ step ].names[ atom ] = names_[ atom ]

            # Copies index of atom to AtomList
            traj[ step ].index[ atom ] = atom
        end

    end

    # Close the poscart.out file
    close( handle_in )

    # Return the vector of AtomList that is the trajectory
    return traj
end
# Reads runtime.inpt, poscart.out, celllens.out and cellangles.out
# - returns a trajectory for the atoms (vector of AtomList) and the cell (vector of Cell_param)
function readPoscarTraj( input_path::T1, poscar_path::T2, cell_length_path::T3, cell_angles_path::T4 ) where { T1 <: AbstractString, T2 <: AbstractString, T3 <: AbstractString, T4 <: AbstractString }
    # Arguments:
    # - input_path: path to the runtime.inpt file (string)
    # - poscar_path: path to the poscart.out file (string)
    # - cell_length_path: path to the celllens.out file (string)
    # - cell_angles_path: path to the cellangles.out file (string)
    # Output
    # - traj: vectors of AtomList - trajectory of atoms
    # - cells: vector of Cell_param - trajectory of the cell

    # Get the species and number of atoms per species
    species, species_nb = pimaim.getSpeciesAndNumber( input_path )

    # Check that it worked
    if species == false || species_nb == false
        # If it fails, returns false, false
        return false, false
    end

    # Read the trajectory of atoms
    traj = readPoscarOut( poscar_path, species, species_nb )
    # Check that it works
    if traj == false
        # If it fails, returns false, false
        return false, false
    end

    # Reads the cell trajectory
    cells = pimaim.readCellParams( cell_length_path, cell_angles_path )
    # Check that it worked
    if cells == false
        # If it did not work, returns false, false
        return false, false
    end

    # Return Trajectory in the shape of
    # - a vector of AtomList for the atoms
    # - a vector of Cell_param for the cells
    return traj, cells
end
# Reads runtime.inpt, poscart.out, cellbox.out
# - returns a trajectory for the atoms (vector of AtomList) and the cell (vector of Cell_param)
function readPoscarTraj( input_path::T1, poscar_path::T2, cell_box_path::T3 ) where { T1 <: AbstractString, T2 <: AbstractString, T3 <: AbstractString }
    # Argument:
    # - input_path: path to runtime.inpt (string)
    # - poscar_path: path to the poscart.out file (string)
    # - cell_box_path: path to cellbox.out (string)
    # Output
    # - traj: vector of AtomList (Trajectory of the atoms)
    # - cells: vector of Cell_param (trajectory of the cell)

    # Get the chemical species and number of element per species
    species, species_nb = pimaim.getSpeciesAndNumber( input_path )

    # Check that all went well
    if species == false || species_nb == false
        # if not, returns false, false
        return false, false
    end

    # Reads the cell trajectory
    cells = pimaim.readCellBox( cell_box_path )

    # Check that cells is ok
    if cells == false
        # If not , returns false, false
        return false, false
    end

    # Reads the atoms trajectory
    traj = readPosCarTraj( poscar_path, species, species_nb )

    # Check that it worked
    if traj == false
        # If not , returns false, false
        return false, false
    end

    # Return Trajectory in the shape of
    # - a vector of AtomList for the atoms
    # - a vector of Cell_param for the cells
    return traj, cells
end
#-----------------------------------------------------------------------------

# Reading positions.out
#-----------------------------------------------------------------------------
# Read positions.out and returns vector of AtomList, requires cell matrices to convert the position to cartesian (converts from Bohr to Angstrom)
function readPositions( path_file::T1, species::Vector{T2}, nb_species::Vector{T3}, cell_matrices::Array{T4,3} ) where { T1 <: AbstractString, T2 <: AbstractString, T3 <: Int, T4 <: Real }
    # Argument:
    # - path_file: path to positions.out (string)
    # - species: names of the species in the cell (vector string)
    # - nb_species: number of elements per chemical species (vector int)
    # - cell_matrices: tensor with the cell trajectory
    # Output
    # - traj: vector of AtomList describing the trajectory
    # Or false if something is wrong with the file

    # Get number of lines
    nb_lines = utils.getNbLines( path_file )

    # If the file does not exists, or is empty, returns false
    if nb_lines == false || nb_lines == 0
        return false
    end

    # Compute number of atoms
    nb_atoms = sum(nb_species)

    # Construct the vector of names of the atoms based on the number of element of species
    # - assumes that atoms are sorted by types
    names_ = atom_mod.buildNames( species, nb_species )

    # Check that the number of lines is coherent with the number of atoms
    if nb_lines % nb_atoms != 0
        # If the format is wrong sends a message and returns false
        print("Problem within file: ",path_file," :\n")
        print("Number of lines is not a multiple of the number of atoms given.\n")
        return false
    end

    # Compute number of steps
    nb_step = Int(nb_lines/nb_atoms)

    # Initialize trajectory
    traj = Vector{AtomList}( undef, nb_step )

    # Open file
    handle_in = open( path_file )

    # Loop over step
    for step=1:nb_step

        # Initialize AtomList for step
        traj[step] = atom_mod.AtomList(nb_atoms)

        # Copy cell matrix for step
        matrix2 = copy( cell_matrices[step,:,:] )

        # Initialize vector for box lengths
        box_len=zeros(Real,3)

        # Compute the lengths
        for i=1:3
            # Compute the norm for all three directions and converts from Bohr to Ang
            box_len[i] = LinearAlgebra.norm( matrix2[:,i] )*conversion.ang2Bohr
        end

        # Realign cell matrix so that its v1 is aligned with z-axis
        cell_matrices[:,:,step] = cell_mod.params2Matrix( cell_mod.matrix2Params( cell_matrices[:,:,step] ) )

        # Loop over atoms
        for atom=1:nb_atoms

            # Initialize vector for positions
            temp = zeros(Real,3)

            # Reads and parse position line for atom
            keys = split( readline( handle_in ) )

            # Loop over dimension
            for i=1:3
                # Converts positions to float and divide by box lengths
                temp[i] = parse(Float64,keys[i])/box_len[i]
            end

            # Converts reduced position to cartesian positions
            # - Loop for dimensions
            for i=1:3
                # - Second loop over dimensions
                for j=1:3
                    # Converts reduced to cartesian
                    traj[step].positions[ i, atom ] += temp[j]*cell_matrices[ i, j, step ]
                end
            end

            # Affects names of atoms to the AtomList
            traj[ step ].names[ atom ] = names_[ atom ]

            # Affects index of atoms to the AtomList
            traj[ step ].index[ atom ] = atom

        end

    end

    # Close atoms
    close( handle_in )

    # Returns trajectory in shape of vector of AtomList
    return traj
end
# Read positions.out and return vector of AtomList, but also works if the format is wrong, requires a cell trajectory in vector of cell_param
function readPositionsCrash( path_file::T1, species::Vector{T2}, nb_species::Vector{T3}, cells::Vector{T4} ) where { T1 <: AbstractString, T2 <: AbstractString, T3 <: Int, T4 <: cell_mod.Cell_param }
    # Argument
    # - path_file: path of positions.out file
    # - species: names of all elements in the structure (vector of string)
    # - nb_species: number of elements for each species (vector of int)
    # - cells: vector of Cell_param for the trajectory of cells
    # Output
    # - traj: vector of AtomList to the trajectory of the atoms
    # OR false if something goes wrong with the file

    # Get the number of lines of the file
    nb_lines = utils.getNbLines( path_file )

    # If the file is empty or does not exists, return false
    if nb_lines == false || nb_lines == 0
        return false
    end

    # Compute number of atoms
    nb_atoms = sum( nb_species )

    # Computes vector of
    names_ = atom_mod.buildNames( species, nb_species )

    # Compute number of steps in the file
    nb_step = Int( trunc( nb_lines/nb_atoms ) )
    # If the number of step is problematic, returns false
    if nb_step < 1
        return false
    end

    # if the size of the cell trajectory is smaller than the number of step
    # reduces the number of steps to that of the number of step in the cell trajectory
    # ( assumption is that the positions didn't write the last steps )
    if size(cells)[1] < nb_step
        nb_step = size(cells)[1]
    end

    # Initialize vector of AtomList for the atom trajectory
    traj = Vector{AtomList}( undef, nb_step )

    # Open the file
    handle_in = open( path_file )

    # Loop over step
    for step=1:nb_step

        # Initialize AtomList for current step
        traj[step] = atom_mod.AtomList(nb_atoms)

        # computes the cell matrix for the current step from a Cell_param
        matrix2 = cell_mod.params2Matrix( cells[step] )

        # Initialize vector for box lengths
        box_len = zeros(Real,3)

        # Loop over dimensions for box length computations and conversion from Bohr to Angstrom
        for i=1:3
            box_len[i] = LinearAlgebra.norm( matrix2[:,i] )*conversion.ang2Bohr
        end

        # Loop over atoms
        for atom=1:nb_atoms

            # Reads and parse the position for current atom
            keys = split( readline( handle_in ) )

            # Initialize temporary file for current position
            temp = zeros(Real,3)

            # Compute position and remove the box length conversion of pimaim
            for i=1:3
                temp[i] = parse(Float64, keys[i] )/box_len[i]
            end

            # Converts positions from reduced to cartesian
            # - Loop over dimension
            for i=1:3
                # - Second Loop over dimensions
                for j=1:3
                    # Converts positions
                    traj[step].positions[ i, atom ] += temp[j]*matrix2[i,j]
                end
            end

            # Affects the names of atoms to the AtomList for current step
            traj[step].names[atom] = names_[atom]

            # Affects the index to the AtomList for current step
            traj[step].index[atom] = atom

        end

    end

    # Close the file
    close( handle_in )

    # Trajectory in the shape of a vector of AtomList
    return traj
end
#-----------------------------------------------------------------------------

#==============================================================================#

# Reading *.xv file
#-----------------------------------------------------------------------------
# Reads structure from .xv file, returns atoms position in AtomList and cell with Cell_param
function readXV( path_file::T1 ) where { T1 <: AbstractString }
    # Argument
    # path_file: path to the *.xv file
    # Output:
    # - atoms: AtomList with atomic informations
    # - cell: Cell_param with information about the cell
    # OR false if something goes wrong

    # Get number of lines
    nb_lines = utils.getNbLines( path_file )
    # If the file does not exists, or is empty, returns false
    if nb_lines == false || nb_lines == 0
        return false
    end

    # Initialize cell matrix
    matrix=zeros(Real,3,3)

    # Open file
    handle_in = open( path_file )

    # Loop over the first three lines (cell matrix)
    for i=1:3
        # Read the line and parse it
        line = split( readline( handle_in ) )
        # Loop over the three first elements
        for j=1:3
            # Converts string to float, and converts lengths from Bohr to Angstroms
            matrix[i,j] = parse(Float64, line[j] )*conversion.bohr2Ang
        end
    end

    # Converts cell matrix to Cell_param
    cell = cell_mod.cellMatrix2Params( matrix )

    # Compute the number of atoms
    nb_atoms = parse(Int64, split( readline( handle_in ) )[1] )

    # Initialize AtomList
    atoms = atom_mod.AtomList( nb_atoms )

    # Loop over atoms
    for atom=1:nb_atoms
        # Read and parse current atom line
        keyword = split( readline( handle_in ) )

        # Get the name of the atom from its atomic number
        atoms.names[ atom ] = periodicTable.z2Names( parse(Int, keyword[2] ) )

        # Use the atomic number as index for AtomList
        atoms.index[ atom ] = atom

        # Loop over the elements 2-5 for the atomic positions
        for i=1:3
            # Parse string to float, and converts position from Bohr to Angstrom
            atoms.positions[ i, atom ] = parse(Float64, keyword[ 2 + i ] )*conversion.bohr2Ang
        end
    end

    # Close file
    close( handle_in )

    # Returns AtomList and Cell_param descrbing the structure
    return atoms, cell
end
#-----------------------------------------------------------------------------

# Reading restart.dat
#-----------------------------------------------------------------------------
# Reads restart.dat file and returns what it contains (positions as AtomList, the rest as vectors or array of real)
# NB: some things are still unclear, may need some work
function readRestart( restart_path::T1, runtime_path::T2 ) where { T1 <: AbstractString, T2 <: AbstractString }
    # Arguments
    # - restart_path: path to the restart.dat file (string)
    # - runtime_path:  path to the runtime.inpt file (string)
    # Output
    # - atoms: AtomList containing positions (optional)
    # - velocity: velocities of the atoms (array, real, nb_atoms,3) (optional)
    # - forces: forces acting on atoms (array, real, nb_atoms, 3) (optional)
    # - polarization: polarization of atoms (vector real)
    # - quadrupolar: barostat parameters (array real)
    # - cell_matrix: cell matrix of the cell
    # - runs_nb_step: vector of int, with number of steps for each precedent
    # OR false if the file does not contain the aforementionned information

    # Get the species and their number of atoms
    species, nb_element_species = getSpeciesAndNumber( runtime_path )

    # Compute total number of different elements
    n_species = size(species)[1]

    # Compute number of atoms
    nb_atoms = sum( nb_element_species )

    # Construct the atom names vector (assume that atoms are sorted by type)
    atoms_names = Vector{AbstractString}(undef, nb_atoms )

    # Loop over species
    for i_spec=1:n_species

        # Compute the offset
        offset_specie = sum( nb_element_species[1:i_spec-1] )

        # Loop over atom species
        for atom_spec=1:nb_element_species[ i_spec ]
            # Compute the species
            atoms_names[ offset_specie + atom_spec ] = species[ i_spec ]
        end
    end

    # Open restart.dat file
    handle_in = open( restart_path )

    # Determine what kind of information is in the file
    # - Whether we have the positions
    check_position  = split( readline( handle_in ) )[1]
    # - Whether we have the velocities
    check_velocity  = split( readline( handle_in ) )[1]
    # - Whether we have the forces
    check_forces    = split( readline( handle_in ) )[1]
    # - Whether we have the polarization information
    check_polarity  = split( readline( handle_in ) )[1]

    # Initialize atoms as false
    atoms = false
    # If there is AtomList if present in file
    if check_position == "T"

        # Initialize AtomList
        atoms = atom_mod.AtomList( nb_atoms )

        # Loop over atoms
        for atom=1:nb_atoms

            # Parse line for current atom
            keyword = split( readline( handle_in ) )

            # Affects index of atom to AtomList
            atoms.index[ atom ] = atom

            # Affects name of atom to AtomList
            atoms.names[ atom ] = atoms_names[ atom ]

            # Loop over dimensions
            for i=1:3
                # Converts from String to Float, converts from bohr to Angstrom
                atoms.positions[ i, atom] = parse(Float64, keyword[i] )*conversion.bohr2Ang
            end
        end
    end

    # Initialize output for velocities
    velocity = false
    # Check if the file contains velocities
    if check_velocity == "T"

        # Initialize vector for velocities
        velocities = zeros(Real, 3, nb_atoms )

        # Loop over the atoms
        for atom=1:nb_atoms

            # Parse the current velocity of atom line
            keyword = split( readline( handle_in ) )

            # Loop over dimension
            for i=1:3
                # converts string to float (may need to convert as well, but unclear about the units of velocities)
                velocities[atom,i] = parse(Float64, keyword[i] )
            end
        end
    end

    # Initialize output for forces
    forces = false

    # Check if the file contains forces
    if check_forces == "T"

        # Initialize forces as Array (real: nb_atoms,3)
        forces = zeros(Real, 3, nb_atoms )

        # Loop over atoms
        for atom=1:nb_atoms

            # Read and parse the force on atom line
            keyword = split( readline( handle_in ) )

            # Loop over dimensions
            for i=1:3
                # Converts strings into floats, may need to convert (but unknown units)
                forces[ i, atom ] = parse(Float64, keyword[i] )
            end
        end
    end

    # Compute the number of previous runs
    nb_runs = parse(Int64, split( readline( handle_in ) )[1])

    # Initialize int vector for number of steps of previous simulations
    runs_nb_step = zeros( nb_runs )

    # Loop over nb_runs following lines to get the number of steps of the previous sim
    for i=1:nb_runs
        # Converts string to int for the number of steps of previous sim
        runs_nb_step[i] = parse(Int64, split( readline(handle_in) )[1] )
    end

    # Initialize polarization vector
    polarization = zeros( 30 )
    # Loop over the next 30 (why 30???) lines
    for i=1:30
        # reads polarization, converts string to float
        polarization[i] = parse(Float64, split( readline(handle_in) )[1] )
    end

    # Initialize quadrupolar matrix (why 2x2??)
    quad = zeros( 2, 2 )
    # Loop over two lines
    for i=1:2

        # Read line and parse it with " " as delimitator
        keys = split( readline( handle_in ) )

        # Loop over two first elements
        for j=1:2
            # Parse from string to float the quadrupolar information
            quad[i,j] = parse(Float64, keys[j] )
        end
    end

    # Ignores three lines (why?)
    for i=1:3
        readline( handle_in )
    end

    # Initialize cell matrix
    cell_matrix = zeros(Real, 3, 3 )
    # Loop over dimension 1
    for i=1:3

        # Reads line and parse with " " as delimitator
        keyword = split( readline( handle_in ) )

        # Loop over dimension 2
        for j=1:3
            # Parse from string to float and send information to cell matrix
            cell_matrix[ i, j ] = parse( Float64, keyword[j] )
        end
    end

    # Initialize box lengths vector
    boxlens = zeros( 3 )

    # Loop over dimensions
    for i=1:3
        # Reads the line and parse with " "
        key = split( readline( handle_in ) )

        # compute box lengths by casting string into float
        boxlens[i] = parse(Float64, key[1] )
    end

    # Converts the cell matrix elements from bohr to angstrom and multiply elements by box length
    # - Loop over dimension 1
    for i=1:3
        # - Loop over dimension 2
        for j=1:3
            # Recompute cell matrix element
            cell_matrix[ i, j ] *= cell_matrix[ i, j ]*boxlens[ i ]*conversion.bohr2Ang
        end
    end

     # returns data
    return atoms, velocity, forces, polarization, quad, cell_matrix, runs_nb_step
end
#------------------------------------------------------------------------------

# Reading Cell informations
#------------------------------------------------------------------------------
# Reads the files celllens.out and cellangles.out and converts them into a vector of Cell_param
function readCellParams( path_file_len::T1, path_file_angles::T2 ) where { T1 <: AbstractString, T2 <: AbstractString }
    # Argument
    # - path_file_len: path to the file celllens.out (string)
    # - path_file_angles: path to the file cellangles.out (string)
    # Output
    # - Vector of Cell_param that gives the information of the trajectory of the cell
    # OR false if something went wrong with the file

    # Get the number of lines in the cellangles file
    nb_lines_angles = utils.getNbLines( path_file_angles )

    # if the file does not exists or is empty, prints a message and returns false
    if nb_lines_angles == false || nb_lines_angles == 0
        print("File cellangles.out does not exists or is empty.\n")
        return false
    end

    # Get the number of lines in the celllens file
    nb_lines_lengths = utils.getNbLines( path_file_len )

    # if the file does not exists or is empty, prints a message and returns false
    if nb_lines_lengths == false
        print("File celllens.out does not exists or is empty.\n")
        return false
    end

    # Check that files have the same number of steps
    nb_step = min( nb_lines_angles, nb_lines_lengths )

    # If not, use the minimum to proceed, but warns the user
    if nb_lines_angles != nb_lines_lengths
        print( "Missmatch between number of lengths and angles, using the minimum to proceed.\n" )
    end

    # Initialize vectors for lenghts and angles
    lengths = zeros(Real, 3, nb_step )
    angles  = zeros(Real, 3, nb_step )

    # Value to convert radiants to degrees
    tau = 180.0/pi

    # Open files
    handle_in_len = open( path_file_len )    # Opens celllens.out
    handle_in_ang = open( path_file_angles ) # Opens cellangles.out

    # Initialize output (vector of Cell_param)
    cells = Vector{ cell_mod.Cell_param }( undef, nb_step )

    # Loop of ver step
    for line=1:nb_step

        # Read  step line ...
        key_len = split( readline( handle_in_len ) ) # for lengths
        key_ang = split( readline( handle_in_ang ) ) # for angles

        # Loop over dimension for conversion
        for i=1:3
            # Converts lengths from Bohr to Angstroms
            lengths[line,i] = parse( Float64, key_len[ 1 + i ] )*conversion.bohr2Ang
            # Converts angles from radians to degrees
            angles[line,i]  = parse( Float64, key_ang[ 1 + i ] )*tau
        end

        # Switch angles around (technicality so that they match the proper order)
        stock              = angles[ line, 1 ]
        angles[ line, 1 ]  = angles[ line, 3 ]
        angles[ line, 3 ]  = stock

        # Puts lengths and angles into the Cell_param for current step
        cells[ line ] = cell_mod.Cell_param( lengths[ line, : ], angles[ line, : ] )
    end

    # Closing files
    close( handle_in_len ) # celllens.out file
    close( handle_in_ang ) # cellangles.out file

    # Returns the vector of Cell_param with cell trajectory
    return cells
end
# Reads the files celllbox.out and converts them into a vector of Cell_param
function readCellBox( path_file::T1 ) where { T1 <: AbstractString }
    # Argument:
    # - path_file: path to the cellbox.out file
    # Output
    # - cells: return a tensor (real; nb_step,3,3) which contains a cell matrix for each step
    # OR false if there is a problem with file

    # Get number of lines of the file
    nb_lines = utils.getNbLines( path_file )

    # If the file is empty or does not exists, prints message and returns false
    if nb_lines == 0 || nb_lines == false
        print("File ",path_file," is empty or does not exist...\n")
        return false
    end

    # Number of line per cell matrix
    nb_line_per_box = 4

    # Check that the format is ok
    if nb_lines % nb_line_per_box != 0
        # If the format is not right, sends a message and returns false
        print("File: ",path_file," does not have the correct number of lines.\n")
        return false
    end

    # Compute number of steps
    nb_step = Int( nb_lines/nb_line_per_box )

    # Initialize cell tensor
    cells = Array{ Real }(undef, 3, 3, nb_step )

    # Opening input file
    handle_in = open( path_file )

    # Loop over steps
    for step=1:nb_step

        # Reading cell reduced matrix for current step
        cells[step,:,:] = zeros(3,3)

        # Loop over dimension
        for i=1:3
            # Read line and parse with " " deliminator
            keys = split( readline( handle_in ) )

            # loop over dimension 2
            for j=1:3
                # parse the element of the matrix from string to float
                cells[step,i,j] = parse( Float64, keys[j] )
            end
        end

        # Reading cell lengths lines and parse with " " as delimitator
        keys = split( readline( handle_in ) )

        # Loop over dimension
        for i=1:3

            # Cast the element from string to float
            length = parse( Float64, keys[i] )

            # Unreduced cell tensor, and converts length from Bohr to Angstrom
            cells[step,:,i] *= length*conversion.bohr2Ang
        end
    end

    # Close file
    close( handle_in )

    # Returns the cell tensor with cell trajectory
    return cells
end
#------------------------------------------------------------------------------

# Reading energy file
#------------------------------------------------------------------------------
# Reading eng1.out file, contains total, potential and kinetic energy
function readEng1( path_file::T1 ) where { T1 <: AbstractString }
    # Argument
    # - path_file: path to the eng1.out file
    # Output
    # - total energy: real vector (nb_step), containing total energy (Kinetic+Potential energy)
    # - kinetic energy: real vector (nb_step), containing kinetic energy
    # - potential energy: real vector (nb_step), containing potential energy
    # OR false, false, false if something is wrong with the file

    # Get number of lines in the input file
    nb_lines = utils.getNbLines( path_file )

    # If the file is empty or does not exists, returns false
    if nb_lines == false || nb_lines == 0
        return false, false, false
    end

    # Initialize output vectors
    enpot = zeros(Real, nb_lines )
    enkin = zeros(Real, nb_lines )

    # Open the input file
    handle_in = open( path_file )

    # Loop over steps
    for step=1:nb_lines
        # Read line, parse with " " as deliminator
        keyword = split( readline( handle_in ) )

        # Converts the string into floats
        enpot[step] = parse(Float64, keyword[2] ) # Second element is the potential energy
        enkin[step] = parse(Float64, keyword[3] ) # Third element is the kinetic energy
    end

    # Close the input file
    close( handle_in )

    # Return the total energy (float), the potential energy (float) and the kinetic energy (float)
    return enpot .+ enkin, enpot, enkin
end
# Reading eng2.out file, contains enthalpy and pv contribution
function readEnthalpy( path_file::T1 ) where { T1 <: AbstractString }
    # Argument:
    # - path_file: path to eng2.out file
    # Output
    # - enthalpy: real vector (nb_step), contains the enthalpy
    # - pv: real vector (nb_step), contains the P*V

    # Get the number of lines in the file
    nb_lines = utils.getNbLines( path_file )

    # If the file does not exists or is empty, returns false
    if nb_lines == false || nb_lines == 0
        return false
    end

    # Initialize output vectors for
    enthalpy = zeros(Real,nb_lines) # - the enthalpy
    pv       = zeros(Real,nb_lines) # - the P*V contribution

    # Open input file
    handle_in = open( path_file )

    # Loop over steps
    for step=1:nb_lines
        # Reading line and parsing using " " as deliminator
        keyword = split(readline(handle_in))

        # Casting strings into float :
        pv[step]       = parse(Float64, keyword[2] )  # - Second column of file is P*V
        enthalpy[step] = parse(Float64, keyword[3] )  # - Third columns is the enthalpy (Total Energy + P*V)
    end

    # Closing file
    close(handle_in)

    # Returning enthalpy and PV contribution
    return enthalpy, pv
end
#------------------------------------------------------------------------------

# Reading thermodynamical data
#------------------------------------------------------------------------------
# Reading pressure.out file, return a real vector contains the pressure
function readPressure( path_file::T1 ) where { T1 <: AbstractString }
    # Argument:
    # - path_file: path to the pressure.out file
    # Output
    # - pressure: vector (real) contains the pressure during trajectory
    # OR false if something is wrong with the file

    # Get number of lines in the file
    nb_lines = utils.getNbLines( path_file )

    # If the file is empty or does not exists, return false
    if nb_lines == false || nb_lines == 0
        return false
    end

    # Initialize vector for output
    pressure = zeros(Real, nb_lines )

    # Open file
    handle_in = open( path_file )

    # Loop over step
    for step=1:nb_lines
        # Read line
        # Parse with " "
        # Get only second element,
        # Cnverts pressure from atomic unit to GPa
        # Converts string to real
        pressure[step] = parse(Float64, split( readline( handle_in ) )[2] )*conversion.au2Gpa
    end

    # Close the file
    close(handle_in)

    # Returns real vector containing pressure
    return pressure
end
# Reading kintemp.out file, returns a real vector containing the temperature
function readTemperature( path_file::T1 ) where { T1 <: AbstractString }
    # Argument:
    # - path_file: path to the kintemp.out file
    # Output
    # - temperature: vector containing the temperature (real; nb_step)
    # OR false if something goes wrong in the file

    # Get number of lines in the file
    nb_lines = utils.getNbLines( path_file )

    # If the file is empty or does not exists, return false
    if nb_lines == 0 || nb_lines == false
        return false
    end

    # Initialize vector for temperature
    temperature = zeros(Real, nb_lines)

    # Open file
    handle_in = open( path_file )

    # Loop over steps
    for step=1:nb_lines
        # Read the line, parse with " ", gets only second element, converts element to float
        temperature[step] = parse(Float64, split( readline(handle_in) )[2] )
    end

    # Close the file
    close( handle_in )

    # Return vector of temperature
    return temperature
end
# Reading cellvol.out returns two vectors containing the volume and some other thing related but unknown
function readVolume( path_file::T1 ) where { T1 <: AbstractString }
    # Argument
    # - path_file: path to the cellvol.out file
    # Output
    # - volume: real vector (nb_step) contains the volume of the cell as a function of time (likely Bohr^3)
    # - other_vol: real vector (nb_step), contains unknown physical value
    # OR false, false if something goes wrong with the file

    # Get the number of lines of the file
    nb_lines = utils.getNbLines( path_file )

    # If the file is empty or does not exists, return false, false
    if nb_lines == 0 || nb_lines == false
        return false, false
    end

    # Initialize vectors for
    volume    = zeros(Real, nb_lines ) # The volume
    other_vol = zeros(Real, nb_lines ) # Whatever the second value in the file is

    # Opens the file
    handle_in = open( path_file )

    # Loop over steps
    for step=1:nb_lines
        # For each value, reads the line, select the element and casts the string into float
        volume[step]    = parse(Float64, split( readline( handle_in ) )[2] ) # For the volume gets the 2nd element
        other_vol[step] = parse(Float64, split( readline( handle_in ) )[3] ) # FOr the other quantity gets the 3rd element
    end

    # Close the file
    close( handle_in )

    # Return the volume and the other quantity
    return volume, other_vol
end
# Reading fulloutput.dat file and returns a data array with thermodynamical values
function readFullOutput( path_file::T1 ) where { T1 <: AbstractString }
    # Argument:
    # - path_file: path to the fulloutput.dat file
    # Output
    # - data: array (nb_step,6) containing thermodynamical values
    # OR false, if something is wrong with the file

    # Gets the number of lines in the file
    nb_lines = utils.getNbLines( path_file )

    # If the file is empty or does not exists, returns false
    if nb_lines == false || nb_lines == 0
        return false
    end

    # Determine data block size
    data_block_size = 7

    # Check the format of the file
    if nb_lines % data_block_size != 0
        # If the format is not right, sends a message to the user and return false
        print("Issue with fulloutput.dat at ", path_file," !\n")
        return false
    end

    # Compute number of steps
    nb_step = Int( nb_lines/7 )

    # Initialize array for output
    data = zeros(nb_step,6)

    # Opens the file
    handle_in = open( path_file )

    # Loop over steps
    for step=1:nb_step

        # Eliminates the first 6 lines that do not contain interesting data
        for i=1:data_block_size-1
            test=readline( handle_in )
        end

        # reads the data line and parse it using " " as delimintator
        keys = split( readline( handle_in ) )

        # Loop over the data elements
        for i=1:6
            # Put data into array (1st column is time)
            data[ step, i ] = parse( Float64, keys[i+1] )
        end
    end

    # Close input file
    close(handle_in)

    # Return the data array
    return data
end
#------------------------------------------------------------------------------

# Extract data from fulloutput.dat
#------------------------------------------------------------------------------
# Extract temperature from a data array from fulloutput.dat
function extractTemperature( data::Array{T1,2} ) where { T1 <: Real }
    # Argument
    # - data: array real (nb_step,6) containing thermodynamical data from fulloutput.dat
    # Output
    # - Real vector (nb_step) containing the temperature (K)

    # Returns the temperature
    return data[:,6]
end
# Extract pressure from a data array from fulloutput.dat
function extractPressure( data::Array{T1,2} ) where { T1 <: Real }
    # Argument
    # - data: array real (nb_step,6) containing thermodynamical data from fulloutput.dat
    # Output
    # - Real vector (nb_step) containing the pressure (GPa)

    # Returns the pressure
    return data[:,5]
end
# Extract volume from a data array from fulloutput.dat
function extractVolume( data::Array{T1,2} ) where { T1 <: Real }
    # Argument
    # - data: array real (nb_step,6) containing thermodynamical data from fulloutput.dat
    # Output
    # - Real vector (nb_step) containing the volume (A^3)

    # Returns the volume
    return data[:,4]
end
# Extract the total energy from a data array from fulloutput.dat
function extractTotalEnergy( data::Array{T1,2} ) where { T1 <: Real }
    # Argument
    # - data: array real (nb_step,6) containing thermodynamical data from fulloutput.dat
    # Output
    # - Real vector (nb_step) containing the total energy (a.u?)

    # Returns the total energy
    return data[:,3]
end
# Extract the kinetic energy from a data array from fulloutput.dat
function extractKineticEnergy( data::Array{T1,2} ) where { T1 <: Real }
    # Argument
    # - data: array real (nb_step,6) containing thermodynamical data from fulloutput.dat
    # Output
    # - Real vector (nb_step) containing the kinetic energy energy (a.u?)

    # Returns the kinetic energy
    return data[:,2]
end
# Extract the potential energy from a data array from fulloutput.dat
function extractPotentialEnergy( data::Array{T1,2} ) where { T1 <: Real }
    # Argument
    # - data: array real (nb_step,6) containing thermodynamical data from fulloutput.dat
    # Output
    # - Real vector (nb_step) containing the potential energy energy (a.u?)

    # Returns the potential energy
    return data[:,1]
end
#------------------------------------------------------------------------------

# Reading f3.dat file containing boroxol fraction as a function of time
#------------------------------------------------------------------------------
function readf3( path_file::T1 ) where { T1 <: AbstractString }
    # Argument:
    # - path_file: path to the f3.dat file
    # Output:
    # - ring_data: contains a real array (nb_step,6) with fraction of various rings
    # OR false if something went wrong with the file

    # Gets the number of lines in the file
    nb_lines = utils.getNbLines( path_file )

    # If the file is empty or does not exists, returns false
    if nb_lines == false || nb_lines == 0
        return false
    end

    # Opens the file
    handle_in = open( path_file )

    # Number of columns in the file
    col_nb = 9

    # Initialize output array
    ring_data = zeros(Real, nb_step, col_nb )

    # skip the first line
    readline( handle_in )

    # Loop over steps
    for step = 1:nb_lines-1
        # Read line for step and parse using " " as deliminator
        keys = split( readline( handle_in ) )

        # Loop over columns
        for col=1:col_nb
            # Converts data from string to real
            ring_data[ step, col ] = parse(Float64, keys[col] )
        end
    end

    # Close the file
    close( handle_in )

    # Return the array with the ring statistics
    return ring_data
end
#------------------------------------------------------------------------------

# Writting data (generic)
#------------------------------------------------------------------------------
# Generic function to write data using an IO handler for the output file
function writeFullData( handle_out::T1, data::Array{T2,2} ) where { T1 <: IO, T2 <: Real }
    # Argument
    # - handle_out: IO handler for the file
    # - data: array containing data (nb_step,n_dim), where nb_step and n_dim are int and left up to the user
    # Output
    # - True, if it works

    # Get the number of steps
    nb_step = size( data )[1]

    # Get the dimension of the data array
    n_dim = size( data )[2]

    # Loop over the array
    for step=1:nb_step
        # Loop over the data array
        for i=1:n_dim
            # Write the data for each dimension
            Base.write( handle_out, string( data[ step, i ], " " ) )
        end

        # Write end of line
        Base.write( handle_out, string( "\n" ) )
    end

    # Return true if it worked
    return true
end
# Generic function to write data using a string for the path of the output file
function writeFullData( file_out::T1, data::Array{T1,2} ) where { T1 <: AbstractString, T2 <: Real }
    # Argument
    # - file_out: path to the output path (string)
    # - data: real array (nb_step,n_dim) containing the data
    # Output
    # - Bool: if the writting was successful or not

    # Open file
    handle_out = open( file_out, "w"  )

    # Writting data
    test = writeFullData( handle_out, data )

    # Check if the data was written
    if ! test
        return false
    end

    # Closing file
    close( handle_out )

    # Return true if it works
    return true
end
#------------------------------------------------------------------------------

# Writing restart.dat file
#------------------------------------------------------------------------------
# Write restart.dat with only positions
function writeRestartPositionsOnly( path_file::T1, atoms::T2, cell::Array{T3,2} ) where { T1 <: AbstractString, T2 <: atom_mod.AtomList, T3 <: Real }
    # Argument
    # - path_file: file to the restart.dat file
    # - atoms: AtomList containing the structure to write
    # - cell: cell_matrix
    # Output
    # - Bool: whether or not the writting went ok

    # Get the number of atoms in the structure
    nb_atoms = atom_mod.getNbAtoms( atoms )

    # Opens output file
    handle_out = open( path_file , "w" )

    # Write the bool concerning the data in the file
    write( handle_out, string( "T\n" ) ) # - We write the position
    write( handle_out, string( "F\n" ) ) # - We don't write velocities
    write( handle_out, string( "F\n" ) ) # - We don't write forces
    write( handle_out, string( "F\n" ) ) # - We don't write polarization

    # Copy the positions from AtomList
    positions_ = copy( atoms.positions )

    # Put the positions into reduced coordinates
    positions_ = cell_mod.getTransformedPosition( positions_, inv(cell.matrix) )

    # Reduced the positions by the cell length (pimaim specific)
    for i=1:3
        positions_[ i,:] *= LinearAlgebra.norm( cell.matrix[ i,:] )
    end

    # Loop over atoms
    for atom=1:nb_atoms
        # Loop over dimensions
        for i=1:3
            # Write positions, converting them into bohr first, and limiting precision to third digit
            write( handle_out, string( round( positions_[ i, atom ]*conversion.ang2Bohr, digits=3 ), " " ) )
        end

        # Writting end of line
        write( handle_out, string("\n") )
    end

    # Copying matrix for manipulation
    matrix = copy( cell )

    # Computing cell lengths
    lengths = cell_mod.cellMatrix2Params(cell).length

    # Loop over dimensions
    for i=1:3
        # Reducing the matrix by the cell matrix length
        matrix[:,i] = matrix[:,i]/lengths[i]
    end

    # Loop over dimension 1
    for i=1:3
        # Loop over dimension 2
        for j=1:3
            # Writting matrix data into file
            write( handle_out, string( round( matrix[i,j], digits=3 ), " ") )
        end

        # Writting end of line
        write( handle_out, string("\n") )
    end

    # Loop over dimensions
    for i=1:3
        # Writting lenghts of the box, converted into bohr into file
        write( handle_out, string( round( LinearAlgebra.norm( cell[:,i])*conversion.ang2Bohr, digits=3 ), "\n" ) )
    end

    # Closing file
    close(handle_out)

    # Returning true if all went ok
    return true
end
#------------------------------------------------------------------------------

# Write Crystal cell
#------------------------------------------------------------------------------
# Writting crystal cell with AtomList and Cell_param
function writeCrystalCell( path_file::T1, atoms::T2, cell::Array{T3,2} ) where { T1 <: AbstractString, T2 <: atom_mod.AtomList, T3 <: Real }
    # Argument:
    # - path_file: path to the restart.dat file to create
    # - atoms: AtomList with positions and other atomic information
    # - cell: cell matrix with information on the cell
    # Output
    # - Bool: whether or not the file was written successfully

    # Get the length of the cell
    lengths = cell_mod.cellMatrix2Params( cell ).length

    # Get the species of element present in the AtomList
    species = atom_mod.getSpecies( atoms )

    # Copying matrix for manipulation
    matrix  = copy( cell )

    # Loop over dimension
    for i=1:3
        # Reducing matrix by the length of cell
        matrix[:,i] = matrix[:,i]/lengths[i]
    end

    # Opening output file
    handle_in = open( path_file, "w" )

    # Loop over dimension 1
    for i=1:3
        # Loop over dimension 2
        for j=1:3
            # Writting matrix element, keeping only precision up to the third digits
            write( handle_in, string( round( matrix[j,i], digits=3 ), " " ) )
        end

        # Write end of line
        write( handle_in, string("\n") )
    end

    # Loop over dimension
    for i=1:3
        # Writing lengths of the cells
        write( handle_in, string( lengths[i], "\n")  )
    end

    # Loop over dimension
    # NB: unsure what this is for
    for i=1:3
        # Writting 1
        write( handle_in, "1\n" )
    end

    # Loop over species
    for i=1:size(species)[1]
        # Writting species Z from name
        write( handle_in, string( periodicTable.names2Z( species[i] ), "\n" ) )

        # Writting arbitrary file name
        write( handle_in, string( species[i], "_quartz.mat\n" ) )
    end

    # Loop over dimension
    for i=1:3
        # Writting length of the cell
        write( handle_in, string( lengths[i], "\n")  )
    end

    # Closing file
    close(handle_in)

    # Returns true if everything works
    return true
end
#------------------------------------------------------------------------------

# Writting acell.txt file
#------------------------------------------------------------------------------
# Writting acell.txt from a cell trajectory using a vector of Cell_param
function writeAcellTxt( path_file::T1, cells::Vector{T2} ) where { T1 <: AbstractString, T2 <: cell_mod.Cell_param }
    # Argument
    # - path_file: path of the output acell.txt file
    # - cells: vector of Cell_param that contains the cell trajectory
    # Output
    # - Bool: whether or not the writting was successful

    # Opening file in writting mode
    handle_out=open( path_file, "w" )

    # Loop over steps
    for step=1:size( cells )[1]
        # Converting current step cell params to cell matrix
        cell_matrix = cell_mod.params2Matrix( cells[ step ] )

        # Loop over dimension 1
        for i=1:3
            # Loop over dimension 2
            for j=1:3
                # Writting cell matrix element
                write( handle_out, string( round(cell_matrix[ i, j ],digits=3), " " ) )
            end
            # Writting end of line
            write( handle_out, "\n" )
        end
    end

    # Close the output file
    close(handle_out)

    # Returns true if all went well
    return true
end
# Writting acell.txt from a cell trajectory using a tensor with cell matrices
function writeAcellTxt( path_file::T1, cells::Array{T2,3} ) where { T1 <: AbstractString, T2 <: Real }
    # Arguments
    # - path_file: path to the acell.txt file to write
    # - cells: tensor (nb_step,3,3) with the cell matrices overtime
    # Output
    # - Bool: whether or not the writting was successful

    # Opening output file
    handle_out=open( path_file, "w" )

    # Loop over steps
    for step = 1:size( cells )[1]
        # Loop over dimensions 1
        for i=1:3
            # Loop over dimensions 2
            for j=1:3
                # Writting cell matrix element, precision up to the third digits
                write( handle_out, string( round( cells[ i, j, step ], digits=3 ), " " ) )
            end
            # Write end of line to file
            write( handle_out, "\n" )
        end
    end

    # Closing the file
    close( handle_out )

    # Return true if the writting was successful
    return true
end
#------------------------------------------------------------------------------

# Writting cellbox.out file
#------------------------------------------------------------------------------
function writeCellbox( path_file::T1, cells::Array{T2,3}, timestep::T3 ) where { T1 <: AbstractString, T2 <: Real, T3 <: Real }
    # Argument
    # - path_file: path to the cellbox.out file to write
    # - cells: tensor describing cell trajectory with cell matrices
    # - timestep: timestep of the simulation (fs)
    # Output
    # - Bool: Whether it successfully wrote the file

    # Opening output file
    handle_out = open( path_file, "w" )

    # Loop over steps
    for step=1:size(cells)[1]
        # Writting step in int and actual simulation time
        write( handle_out, string( step, " ", step*timestep, " " ) )

        # Loop over dimension 1
        for i=1:3
            # Loop over dimension 2
            for j=1:3
                # Write cell matrix element
                write( handle_out, string( cells[ i, j, step ], " " ) )
            end
        end

        # Write the total volume of the cell
        write( handle_out, string( " ", det( cells[ :, :, step ] ), "\n" ) )
    end

    # Closing file
    close(handle_out)

    # Returns true if it worked
    return true
end
#------------------------------------------------------------------------------

end
