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
                atoms.positions[atom,i] += matrix2[j,i]*temp[j]
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
                traj[step].positions[atom,i] = parse( Float64, keys[i] )*conversion.bohr2Ang
            end

            # Put the name of the atoms in
            traj[step].names[atom] = names_[atom]

            # Put the index of the atoms in
            traj[step].index[atom] = atom
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
    nb_step = Int(nb_lines/nb_atoms)

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
                traj[step].positions[atom,i] = parse( Float64, keys[i] )*conversion.bohr2Ang
            end

            # conversion of positions
            positions_temp = matrix2 * inv_matrix * traj[step].positions[atom,:]

            # Copies atomic positions to AtomList
            traj[step].positions[atom,:] = positions_temp

            # Copies names of atom to AtomList
            traj[step].names[atom] = names_[atom]

            # Copies index of atom to AtomList
            traj[step].index[atom] = atom
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
function readPositions( path_file::T1, species::Vector{T2}, nb_species::Vector{T3}, cell_matrices::Array{T4,3} ) where { T1 <: AbstractString, T2 <: AbstractString, T3 <: Int, T4 <: Real }
    # Argument:
    # - path_file:
    # - species:
    # - nb_species:
    # - cell_matrices
    # Output
    # - traj: vector of AtomList describing the trajectory

    nb_lines = utils.getNbLines( path_file )
    if nb_lines == false
        # utils.getNbLines will have return the error message
        return false
    end

    nb_atoms = sum(nb_species)
    names_=atom_mod.buildNames( species, nb_species )
    if nb_lines % nb_atoms != 0
        print("Problem within file: ",path_file," :\n")
        print("Number of lines is not a multiple of the number of atoms given.\n")
        return false
    end

    nb_step = Int(nb_lines/nb_atoms)

    traj = Vector{AtomList}( undef, nb_step )

    handle_in = open( path_file )

    for step=1:nb_step

        box_len=zeros(Real,3)

        traj[step] = atom_mod.AtomList(nb_atoms)

        matrix2 = copy( cell_matrix[step,:,:] )

        for i=1:3
            box_len[i] = LinearAlgebra.norm( matrix2[:,i] )*conversion.ang2Bohr
        end

        cells[step].matrix = cell_mod.params2Matrix( cell_mod.cellMatrix2Params( cells[step] ) ).matrix

        for atom=1:nb_atoms

            temp = zeros(Real,3)

            keys = split( readline( handle_in ) )

            for i=1:3
                temp[i] = parse(Float64,keys[i])/box_len[i]
            end

            for i=1:3
                for j=1:3
                    traj[step].positions[atom,i] += temp[j]*cells[step].matrix[i,j]
                end
            end

            traj[step].names[atom] = names_[atom]
            traj[step].index[atom] = atom

        end

    end

    close( handle_in )

    return traj
end
function readPositionsCrash( path_file::T1, species::Vector{T2}, nb_species::Vector{T3}, cells::Vector{T4} ) where { T1 <: AbstractString, T2 <: AbstractString, T3 <: Int, T4 <: cell_mod.Cell_param }

    #---------------------------------------------------------------------
    nb_lines = utils.getNbLines( path_file )
    if nb_lines == false
        return false
    end
    #---------------------------------------------------------------------

    #---------------------------------------------------------------------
    nb_atoms = sum(nb_species)
    names_=atom_mod.buildNames( species, nb_species )
    nb_step = Int( trunc( nb_lines/nb_atoms ) )-1
    if nb_step < 1
        return false
    end
    if size(cells)[1] < nb_step
        nb_step = size(cells)[1]
    end
    traj = Vector{AtomList}( undef, nb_step )
    #---------------------------------------------------------------------

    #---------------------------------------------------------------------
    handle_in = open( path_file )
    for step=1:nb_step
        box_len=zeros(Real,3)
        traj[step] = atom_mod.AtomList(nb_atoms)
        matrix2 = copy( cell_mod.params2Matrix( cells[step] ).matrix )
        for i=1:3
            box_len[i] = LinearAlgebra.norm( matrix2[:,i] )*conversion.ang2Bohr
        end
        for atom=1:nb_atoms
            temp = zeros(Real,3)
            keys = split( readline( handle_in ) )
            for i=1:3
                temp[i] = parse(Float64,keys[i])/box_len[i]
            end
            for i=1:3
                for j=1:3
                    traj[step].positions[atom,i] += temp[j]*matrix2[i,j]
                end
            end
            traj[step].names[atom] = names_[atom]
            traj[step].index[atom] = atom
        end
    end
    #---------------------------------------------------------------------
    close( handle_in )

    return traj
end
#-----------------------------------------------------------------------------

#==============================================================================#

# *.xv
#-----------------------------------------------------------------------------
function readXV( path_file::T1 ) where { T1 <: AbstractString }

    matrix=zeros(Real,3,3)

    handle_in = open( path_file )
    for i=1:3
        line = split(readline( handle_in ))
        for j=1:3
            matrix[i,j] = parse(Float64, line[j] )*conversion.bohr2Ang
        end
    end

    cell = cell_mod.cellMatrix2Params( matrix )

    nb_atoms = parse(Int64, split( readline( handle_in ) )[1] )

    atoms = atom_mod.AtomList( nb_atoms )

    for atom=1:nb_atoms
        keyword = split( readline( handle_in ) )
        atoms.names[atom] = periodicTable.z2Names( parse(Int, keyword[2] ) )
        atoms.index[atom] = atom
        for i=1:3
            atoms.positions[atom,i] = parse(Float64, keyword[2+i] )*conversion.bohr2Ang
        end
    end
    close( handle_in )
    return atoms, cell
end
#-----------------------------------------------------------------------------

#  restart.dat
#-----------------------------------------------------------------------------
function readRestart( restart_path::T1, runtime_path::T2 ) where { T1 <: AbstractString, T2 <: AbstractString }

    #-----------------------------------------------------------------------
    species, nb_element_species = getSpeciesAndNumber( runtime_path )
    n_species = size(species)[1]
    nb_atoms = sum( nb_element_species )
    atoms_names = Vector{AbstractString}(undef,nb_atoms)
    for i_spec=1:n_species
        offset_specie = sum( nb_element_species[1:i_spec-1] )
        for atom_spec=1:nb_element_species[i_spec]
            atoms_names[ offset_specie + atom_spec ] = species[ i_spec ]
        end
    end
    #-----------------------------------------------------------------------

    #-----------------------------------------------------------------------
    handle_in = open( restart_path )
    check_position  = split( readline( handle_in ) )[1]
    check_velocity  = split( readline( handle_in ) )[1]
    check_forces    = split( readline( handle_in ) )[1]
    check_polarity  = split( readline( handle_in ) )[1]
    #-----------------------------------------------------------------------

    #-----------------------------------------------------------------------
    atoms = false
    if check_position == "T"
        atoms = atom_mod.AtomList( nb_atoms )
        for atom=1:nb_atoms
            keyword = split( readline( handle_in ) )
            atoms.index[atom] = atom
            atoms.names[atom] = atoms_names[atom]
            for i=1:3
                atoms.positions[atom,i] = parse(Float64, keyword[i] )*conversion.bohr2Ang
            end
        end
    end
    #-----------------------------------------------------------------------

    #-----------------------------------------------------------------------
    velocity = false
    if check_velocity == "T"
        velocities = zeros(Real, nb_atoms, 3)
        for atom=1:nb_atoms
            keyword = split( readline( handle_in ) )
            for i=1:3
                velocities[atom,i] = parse(Float64, keyword[i] )
            end
        end
    end
    #-----------------------------------------------------------------------

    #-----------------------------------------------------------------------
    forces = false
    if check_forces == "T"
        forces = zeros(Real, nb_atoms, 3)
        for atom=1:nb_atoms
            keyword = split( readline( handle_in ) )
            for i=1:3
                forces[ atom, i ] = parse(Float64, keyword[i] )
            end
        end
    end
    #-----------------------------------------------------------------------

    #-----------------------------------------------------------------------
    nb_runs = parse(Int64, split( readline( handle_in ) )[1])
    runs_nb_step = zeros( nb_runs )
    for i=1:nb_runs
        runs_nb_step[i] = parse(Int64, split( readline(handle_in) )[1] )
    end
    #-----------------------------------------------------------------------

    #-----------------------------------------------------------------------
    polarization = zeros( 30 )
    for i=1:30
        polarization[i] = parse(Float64, split( readline(handle_in) )[1] )
    end
    #-----------------------------------------------------------------------

    #-----------------------------------------------------------------------
    quad = zeros( 2, 2 )
    for i=1:2
        keys = split( readline( handle_in ) )
        for j=1:2
            quad[i,j] = parse(Float64, keys[j] )
        end
    end
    #-----------------------------------------------------------------------

    #-----------------------------------------------------------------------
    for i=1:3
        readline( handle_in )
    end
    #-----------------------------------------------------------------------

    #-----------------------------------------------------------------------
    cell_matrix = zeros(Real,3,3)
    for i=1:3
        keyword = split( readline(handle_in) )
        for j=1:3
            cell_matrix[i,j] = parse( Float64, keyword[j] )
        end
    end
    boxlens = zeros(3)
    for i=1:3
        key = split( readline( handle_in ) )
        boxlens[i] = parse(Float64, key[1] )
    end
    for i=1:3
        for j=1:3
            cell_matrix[i,j] *= cell_matrix[i,j]*boxlens[i]*conversion.bohr2Ang
        end
    end
    #-----------------------------------------------------------------------

    return atoms, velocity, forces, polarization, quad, cell_matrix, runs_nb_step
end
#------------------------------------------------------------------------------

# Reading Cell informations
#------------------------------------------------------------------------------
# celllens.out, cellangles.out -> Cell_param
function readCellParams( path_file_len::T1, path_file_angles::T2 ) where { T1 <: AbstractString, T2 <: AbstractString }

    #-------------------------------------------------
    nb_lines_angles = utils.getNbLines( path_file_angles )
    if nb_lines_angles == false
        return false
    end
    nb_lines_lengths = utils.getNbLines( path_file_len )
    if nb_lines_lengths == false
        return false
    end
    nb_step=min( nb_lines_angles, nb_lines_lengths )
    if nb_lines_angles != nb_lines_lengths
        print("Missmatch between number of lengths and angles, using the minimum to proceed.\n")
    end
    #-------------------------------------------------

    #-------------------------------------------------
    lengths = zeros(Real, nb_step , 3 )
    angles  = zeros(Real, nb_step , 3 )
    #-------------------------------------------------

    #-------------------------------------------------
    tau=180/pi
    handle_in_len = open( path_file_len )
    handle_in_ang = open( path_file_angles )
    cells = Vector{ cell_mod.Cell_param }( undef, nb_step )
    for line=1:nb_step
        key_len = split( readline( handle_in_len ) )
        key_ang = split( readline( handle_in_ang ) )
        for i=1:3
            lengths[line,i] = parse( Float64, key_len[1+i] )*conversion.bohr2Ang
            angles[line,i]  = parse( Float64, key_ang[1+i] )*tau
        end
        stock=angles[line,1]
        angles[line,1]=angles[line,3]
        angles[line,3]=stock
        cells[ line ] = cell_mod.Cell_param( lengths[line,:], angles[line,:] )
        #angles[line,:]=circshift(angles[line,:],1)
    end
    close( handle_in_len )
    close( handle_in_ang )
    #-------------------------------------------------

    return cells
end
# cellbox.out -> Cell_matrix
function readCellBox( path_file::T1 ) where { T1 <: AbstractString }
    nb_lines = utils.getNbLines( path_file )

    if nb_lines == 0
        print("File ",path_file," is empty...\n")
        return false
    end

    nb_line_per_box = 4
    if nb_lines == false
        return false
    end
    if nb_lines % nb_line_per_box != 0
        print("File: ",path_file," does not have the correct number of lines.\n")
        return false
    end
    nb_step = Int( nb_lines/nb_line_per_box )

    cells = Array{ Real }(undef, nb_step, 3, 3 )
    handle_in = open( path_file )
    for step = 1:nb_step
        # Reading cell reduced matrix
        cells[step,:,:] = zeros(3,3)
        for i=1:3
            keys = split( readline( handle_in ) )
            for j=1:3
                cells[step,i,j] = parse( Float64, keys[j] )
            end
        end
        # Reading cell lengths
        keys = split( readline( handle_in ) )
        for i=1:3
            length = parse( Float64, keys[i] )
            cells[step,:,i] *= length*conversion.bohr2Ang
        end
    end
    close( handle_in )

    return cells
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
function readEnergy( path_file::T1 ) where { T1 <: AbstractString }

    #---------------------------------------------------------------------
    nb_lines = utils.getNbLines( path_file )
    if nb_lines == false
        return false
    end
    #---------------------------------------------------------------------

    enpot = zeros(Real,nb_lines)
    enkin = zeros(Real,nb_lines)
    handle_in = open( path_file )
    for step=1:nb_lines
        keyword = split(readline(handle_in))
        enpot[step] = parse(Float64, keyword[2] )
        enkin[step] = parse(Float64, keyword[3] )
    end
    close(handle_in)

    return enpot.+enkin, enpot, enkin
end
function readEnthalpy( path_file::T1 ) where { T1 <: AbstractString }

    #---------------------------------------------------------------------
    nb_lines = utils.getNbLines( path_file )
    if nb_lines == false
        return false
    end
    #---------------------------------------------------------------------

    enthalpy = zeros(Real,nb_lines)
    pv = zeros(Real,nb_lines)
    handle_in = open( path_file )
    for step=1:nb_lines
        keyword = split(readline(handle_in))
        pv[step] = parse(Float64, keyword[2] )
        enthalpy[step] = parse(Float64, keyword[3] )
    end
    close(handle_in)

    return enthalpy, pv
end
function readPressure( path_file::T1 ) where { T1 <: AbstractString }

    #---------------------------------------------------------------------
    nb_lines = utils.getNbLines( path_file )
    if nb_lines == false
        return false
    end
    #---------------------------------------------------------------------

    pressure = zeros(Real,nb_lines)
    handle_in = open( path_file )
    for step=1:nb_lines
        pressure[step] = parse(Float64, split(readline(handle_in))[2] )*conversion.au2Gpa
    end
    close(handle_in)

    return pressure
end
function readTemperature( path_file::T1 ) where { T1 <: AbstractString }
    nb_lines = utils.getNbLines( path_file )
    if nb_lines == 0 || nb_lines == false
        return false
    end

    temperature = zeros(Real, nb_lines)
    handle_in = open( path_file )
    for step=1:nb_lines
        temperature[step] = parse(Float64, split( readline(handle_in) )[2] )
    end
    close( handle_in )
    return temperature
end
function readVolume( path_file::T1 ) where { T1 <: AbstractString }

    nb_lines = utils.getNbLines( path_file )
    if nb_lines == 0 || nb_lines == false
        return false
    end

    volume = zeros(Real, nb_lines)
    handle_in = open( path_file )
    for step=1:nb_lines
        volume[step] = parse(Float64, split( readline(handle_in) )[2] )
    end
    close( handle_in )
    return volume
end
# fulloutput.dat
function readFullOutput( path_file::T1 ) where { T1 <: AbstractString }
    nb_lines = utils.getNbLines( path_file)
    if nb_lines == false
        return false
    end
    if nb_lines % 7 != 0
        print("Issue with fulloutput.dat at ", path_file,"!\n")
        return false
    end
    nb_step = Int( nb_lines/7 )
    handle_in = open( path_file )
    data = zeros(nb_step,6)
    for step=1:nb_step
        for i=1:6
            test=readline( handle_in )
        end
        keys = readline( handle_in )
        for i=1:6
            data[ step, i ] = parse( Float64, split(keys)[i+1] )
        end
    end
    close(handle_in)
    return data
end
#------------------------------------------------------------------------------

# f3.dat
#------------------------------------------------------------------------------
function readf3( path_file::T1 ) where { T1 <: AbstractString }
    handle_in = open( path_file )
    nb_step = 0
    while ! eof( handle_in )
        readline( handle_in )
        nb_step += 1
    end
    seekstart( handle_in )
    col_nb = 9
    readline( handle_in )
    nb_step = nb_step-1
    ring_data = zeros( nb_step, col_nb )
    for step = 1:nb_step
        keys = split( readline( handle_in ) )
        for col=1:col_nb
            ring_data[ step, col ] = parse(Float64, keys[col] )
        end
    end
    close( handle_in )
    return ring_data
end
#------------------------------------------------------------------------------

# Extract from fulloutput.dat data
#------------------------------------------------------------------------------
function extractTemperature( data::Array{T1,2} ) where { T1 <: Real }
    return data[:,6]
end
function extractPressure( data::Array{T1,2} ) where { T1 <: Real }
    return data[:,5]
end
function extractVolume( data::Array{T1,2} ) where { T1 <: Real }
    return data[:,4]
end
function extractTotalEnergy( data::Array{T1,2} ) where { T1 <: Real }
    return data[:,3]
end
function extractKineticEnergy( data::Array{T1,2} ) where { T1 <: Real }
    return data[:,2]
end
function extractPotentialEnergy( data::Array{T1,2} ) where { T1 <: Real }
    return data[:,1]
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
function writeFullData( handle_out::T1, data::Array{T2,2} ) where { T1 <: IO, T2 <: Real }
    nb_step = size( data )[1]
    n_dim = size( data )[2]
    for step=1:nb_step
        for i=1:n_dim
            Base.write( handle_out, string( data[step,i], " " ) )
        end
        Base.write( handle_out, string( "\n" ) )
    end
    return true
end
function writeFullData( file_out::T1, data::Array{T1,2} ) where { T1 <: AbstractString, T2 <: Real }
    handle_out = open( file_out, "w"  )
    test=writeFullData( handle_out, data )
    close( handle_out )
    return test
end
function writeRestart( path_file::T1, atoms::T2, cell::Array{T3,2} ) where { T1 <: AbstractString, T2 <: atom_mod.AtomList, T3 <: Real }
    nb_atoms=atom_mod.getNbAtoms(atoms)
    handle_out=open( path_file , "w" )
    write( handle_out, string("T\n") )
    write( handle_out, string("F\n") )
    write( handle_out, string("F\n") )
    write( handle_out, string("F\n") )
    positions_ = copy(atoms.positions)
    positions_ = cell_mod.getTransformedPosition( positions_, inv(cell.matrix) )
    for i=1:3
        positions_[:,i] *= LinearAlgebra.norm( cell.matrix[:,i] )
    end
    for atom=1:nb_atoms
        for i=1:3
            write( handle_out, string( round(positions_[atom,i]*conversion.ang2Bohr,digits=3)," " ) )
        end
        write( handle_out, string("\n") )
    end

    matrix = copy( cell )

    lengths = cell_mod.cellMatrix2Params(cell).length
    for i=1:3
        matrix[:,i] = matrix[:,i]/lengths[i]
    end

    for i=1:3
        for j=1:3
            write( handle_out, string( round(matrix[i,j],digits=3), " ") )
        end
        write( handle_out, string("\n") )
    end

    for i=1:3
        write(handle_out,string(round(LinearAlgebra.norm(cell[:,i])*conversion.ang2Bohr,digits=3),"\n"))
    end

    close(handle_out)

    return true
end
function writeCrystalCell( path_file::T1, atoms::T2, cell::Array{T3,2} ) where { T1 <: AbstractString, T2 <: atom_mod.AtomList, T3 <: Real }

    lengths = cell_mod.cellMatrix2Params(cell).length

    species = atom_mod.getSpecies( atoms )

    matrix  = copy( cell )

    for i=1:3
        matrix[:,i] = matrix[:,i]/lengths[i]
    end

    handle_in = open( path_file, "w" )

    for i=1:3
        for j=1:3
            write( handle_in, string( round(matrix[j,i],digits=3), " " ) )
        end
        write( handle_in, string("\n") )
    end

    for i=1:3
        write( handle_in, string( lengths[i], "\n")  )
    end

    for i=1:3
        write( handle_in, "1\n" )
    end

    for i=1:size(species)[1]
        write( handle_in, string( periodicTable.names2Z( species[i] ), "\n" ) )
        write( handle_in, string( species[i], "_quartz.mat\n" ) )
    end

    for i=1:3
        write( handle_in, string( lengths[i], "\n")  )
    end

    close(handle_in)

    return true
end
function writeAcellTxt( path_file::T1, cells::Vector{T2} ) where { T1 <: AbstractString, T2 <: cell_mod.Cell_param }

    handle_out=open( path_file, "w" )

    for step=1:size(cells)[1]

        cell_matrix = cell_mod.params2Matrix(cells[step])
        for i=1:3
            for j=1:3
                write( handle_out, string( round(cell_matrix[i,j],digits=3), " " ) )
            end
            write( handle_out, "\n" )
        end
    end

    close(handle_out)

    return true
end
function writeAcellTxt( path_file::T1, cells::Array{T2,3} ) where { T1 <: AbstractString, T2 <: Real }


    handle_out=open( path_file, "w" )

    for step = 1:size(cells)[1]
        for i=1:3
            for j=1:3
                write( handle_out, string( round(cells[step,i,j],digits=3), " " ) )
            end
            write( handle_out, "\n" )
        end
    end

    close(handle_out)

    return true
end
function writeCell( path_file::T1, cells::Array{T2,3}, timestep::T3 ) where { T1 <: AbstractString, T2 <: Real, T3 <: Real }

    handle_out = open( path_file, "w" )

    for step=1:size(cells)[1]

        write( handle_out, string( step, " ", step*timestep, " " ) )

        for i=1:3
            for j=1:3

                write( handle_out, string( cells[step,i,j], " " ) )
            end
        end

        write( handle_out, string( " ", det( cells[step,:,:] ), "\n" ) )
    end

    close(handle_out)

    return true
end
#------------------------------------------------------------------------------

end
