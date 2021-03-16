module pdb

# Loading necessary modules from LibAtomSim
using utils
using atom_mod
using cell_mod
using periodicTable

# Exporting modules functions
export getNbSteps, getNbAtoms
export readStructure, readTraj, readStructureAtomList, readStructureAtomMolList, readTrajAtomListFixedCell, readTrajAtomList
export addJustifyRight, addJustifyLeft
export writeCRYST1, writeMODEL
export writePdb, writePdbPivClustering, writePdbPlumed

# Descriptors
# - Set of functions to deal with pdb files

# Get number of atoms and steps
#-------------------------------------------------------------------------------
# Get number of steps of the *.pdb file
function getNbSteps( file_path::T1 ) where { T1 <: AbstractString }
    # Argument
    # - file_path: path of the *.pdb file
    # Output
    # - nb_step: number of step
    # OR false if something went wrong with file

    # Check that the file exists
    if ! isfile( file_path )
        # If it does not, sends a message and returns false
        print( "No pdb file found at ", file_path," !\n" )
        return false
    end

    # Initialize step counter (two counter to ensure against file corruption)
    nb_step  = 0
    nb_step2 = 0

    # Opens input file
    handle_in = open( file_path )

    # Loop as long as possible over file lines
    while !eof( handle_in )
        # Read current line and parse with " " deliminator
        keyword1 = split(readline(handle_in))[1]

        # first counter is triggered by the END keyword
        if keyword1 == "END"
            nb_step += 1

        # second counter is triggered by "CRYST1" keyword
        elseif keyword1 == "CRYST1"
            nb_step2 +=1
        end
    end

    # Close input file
    close( handle_in )

    # Check that both counter exists (file is not corrupted)
    if nb_step == nb_step2
        # Returns number of steps
        return nb_step
    else
        # If problem, sends a message and returns false
        print("PDB file at ",file_path," is probably corrupted.\n")
        return false
    end
end
# Get number of atoms
function getNbAtoms( file_path::T1 ) where { T1 <: AbstractString }
    # Argument
    # - file_path: path to the input file (string)
    # Output
    # - nb_atoms: number of atoms in the structure (int)
    # OR false, if something went wrong


    # Check that file exists
    if ! isfile( file_path )
        # If not, send an error message, and returns false
        print("No pdb file found at ",file_path," !\n")
        return false
    end

    # Initialize number of atoms for output
    nb_atoms  = 0

    # Open input file
    handle_in = open( file_path )

    # Reading in loop
    while ! eof( handle_in )
        # Reads and parse line, selecting only first column of file
        keyword1 = split( readline( handle_in ) )[1]

        # Increment counter as long as we encounter the "ATOM" keyword
        if keyword1 == "ATOM"
            nb_atoms += 1
        # Stops reading when reading "END" keyword
        elseif keyword1 == "END"
            break
        end
    end

    # Close input file
    close( handle_in )

    # Returns number of atoms
    return nb_atoms
end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Reads a .pdb file containing a single structure
function readStructure( file_path::T1 ) where { T1 <: AbstractString }
    # Argument
    # - file_path: path of the pdb file
    # Output
    # - names: vector (nb_atoms) of string with atomic names
    # - positions: array (3,nb_atoms), contains the atomic positions
    # - cell: Cell_param that contains all cell informations
    # OR false, false, false if something went wrong

    # Get number of atoms and check that file exists
    nb_atoms = getNbAtoms(file_path)

    # If file does not exists returns false, false
    if nb_atoms == false
        return false, false, false
    end

    # Opens input file
    handle_in = open( file_path )

    # Initialize vectors for lengths and angles of cell
    lengths = zeros(Real, 3 )
    angles  = zeros(Real, 3 )

    # Reads and parse first line
    keyword = split( readline( handle_in ) )

    # Reading cell information from CRYST1 line
    if keyword[1] == "CRYST1"
        # Loop over dimensions
        for i=1:3
            # Lengths are elements 2-4
            lengths[i] = parse( Float64, keyword[i+1] )

            # Angles are elements 5-7
            angles[i] = parse( Float64, keyword[i+4] )
        end

    # If initial element is not CRYST1, return error
    else
        print( "Problem with .pdb file at: ", file_path, "\n" )
        return false, false, false
    end

    # Compact all information into a Cell Param object
    cell = cell_mod.Cell_param( lengths, angles )

    # Initialize position arrays
    positions = zeros(Real, 3, nb_atoms )

    # Initialize vector for atom names
    names = Vector{AbstractString}(undef, nb_atoms )

    # Loop over atoms
    for atom=1:nb_atoms
        # Reading line
        keyword = split( readline( handle_in ) )

        # Checking that the first key is ATOM
        if keyword[1] == "ATOM"
            # Atom name is the third element of the line
            names[atom] = keyword[3]

            # Loop over dimensions
            for i=1:3
                # Positions are elements 5-7 in the ATOM lines
                positions[ i, atom ] = parse(Float64, keyword[ 5 + i ] )
            end

        # If not, we have a problem and return false, false
        else
            print( "Problem reading .pdb file at: ", file_path, " at ATOM keyword.\n" )
            return false, false, false
        end
    end

    # Closes input file
    close( handle_in )

    # We return the names of atoms, their position and the cell information
    return names, positions, cell
end
# Reads a .pdb file for a trajectory
function readTraj( file_path::T1 ) where { T1 <: AbstractString }
    # Argument
    # - file_path: path to the input file (string)
    # Output
    # - names: vector of string with atoms names
    # - positions: array (3,nb_atoms,nb_step) with atomic positions
    # - cells: vector of cell param with cell information
    # OR false, false, false if something went wrong

    # Getting number of steps
    nb_step = getNbSteps( file_path )

    if nb_step == false
        return false, false, false
    end

    # Get number of atoms
    nb_atoms = getNbAtoms( file_path )

    # Check that we could get an number of atom
    if nb_atoms == false
        return false, false, false
    end

    # Open input file
    handle_in = open( file_path )

    # Initialize vector for atom names
    names = Vector{AbstractString}(undef, nb_atoms )

    # Initialize array for atomic positions
    positions = zeros(Real, 3, nb_atoms, nb_step )

    # Initialize cell param vector
    cells = Vector{cell_mod.Cell_param}(undef,nb_step)

    # Loop over steps
    for step=1:nb_step
        # If "CRYST1" keyword, reads cell information
        if split(readline( handle_in ) )[1] == "CRYST1"
            # Initialize vectors for cell lengths and angles for current step
            lengths = zeros(Real, 3 )
            angles  = zeros(Real, 3 )

            # Loop over dimensions
            for i=1:3
                # Lengths of cells are columns 2-4 (+ parse to float)
                lengths[i] = parse( Float64, keyword[ i + 1 ] )

                # Angles of cells are columns 5-7 (+ parse to float)
                angles[i]  = parse( Float64, keyword[ i + 4 ] )
            end

            # Assemble cell lengths and angles into CellParam
            cells[step] = cell_mod.Cell_param( lengths, angles)
        # If not returns error message and false, false, false
        else
            print("Problem reading file ", file_path, " at step: ", step, "\n" )
            return false, false, false
        end

        # Loop over atoms
        for atom=1:nb_atoms
            # Read "ATOM" line
            keyword = split( readline( handle_in ) )

            # If "ATOM" keyword
            if keyword[1] == "ATOM"
                # At first step only...
                if step == 1
                    # Reads the atom names as third column
                    names[ count_ ] = keyword[3]
                end

                # Loop over dimensions
                for i=1:3
                    # Read dimensions as columns 6-8
                    positions[ i, count_ ] = parse( Float64, keyword[ 5 + i ] )
                end
            end
        end

        # Reads "END" line
        readline( handle_in )
    end

    # Returns names, positions and cell informations
    return names, positions, cells
end
# Reads a .pdb file containing a single structure
function readStructureAtomList( file_path::T1 ) where { T1 <: AbstractString }
    # Argument
    # - file_path: path to the input file
    # Output
    # - atoms: AtomList, contains atomic information
    # - cells: Cell_param, contains cell information
    # OR false, false

    # Get number of atoms
    nb_atoms = getNbAtoms(file_path)

    # If something went wrong, returns false, false (message was sent by getNbAtoms)
    if nb_atoms == false
        return false, false
    end

    # Open input file
    handle_in = open( file_path )

    # Initialize vectors for cell lengths and angles
    lengths = zeros( Real, 3 )
    angles  = zeros( Real, 3 )

    # Read first line and parse it with " " deliminator
    keyword=split( readline( handle_in ) )

    # Check that the first element is "CRYST1" (cell information)
    if keyword[1] == "CRYST1"
        # Loop over dimensions
        for i=1:3
            # Casts elements 2-4 into float as lengths of cell
            lengths[i] = parse( Float64, keyword[ i + 1 ] )

            # Casts elements 5-7 into float as angles of cell
            angles[i]  = parse( Float64, keyword[ i + 4 ] )
        end

    # If there is no "CRYST1" keyword, there is a problem, returns false, false
    else
        print("Problem with .pdb file at: ",file_path,"\n")
        return false, false
    end

    # Assemble cell information into Cell_param
    cell = cell_mod.Cell_param( lengths, angles )

    # Initialize AtomList for atomic informations
    atoms = atom_mod.AtomList( nb_atoms )

    # Loop over atoms
    for atom=1:nb_atoms
        # Read line and parse with " " deliminator
        keyword = split( readline( handle_in ) )

        # Check that first column starts with "ATOM"
        if keyword[1] == "ATOM"
            # Parse third column as element name
            atoms.names[atom] = keyword[3]

            # Get index of atom
            atoms.index[atom]  = atom

            # Loop over dimensions
            for i=1:3
                # Parse elements 6-8 as particle positions
                atoms.positions[ i, atom ] = parse(Float64, keyword[ 5 + i ] )
            end

        # If column does not start with "ATOM", there is a problem with the file
        else
            print("Problem with .pdb file at: ",file_path," at ATOM keyword.\n")
            return false
        end
    end

    # Close input file
    close( handle_in )

    # Returns atomic (atoms) and cell information (cell)
    return atoms, cell
end
# Reads a .pdb file containing a single structure, returns AtomMolList
function readStructureAtomMolList( file_path::T1 ) where { T1 <: AbstractString }
    # Argument
    # - file_path: path to the input file
    # Output
    # - atoms: AtomMolList for atomic information
    # - cell: Cell_param for cell information

    # Get number of atoms
    nb_atoms = getNbAtoms(file_path)

    # If something went wrong, returns false, false
    if nb_atoms == false
        return false, false
    end

    # Opens input file
    handle_in = open( file_path )

    # Initialize Cell_param for cell info
    cell = cell_mod.Cell_param()

    # Read first line and parse with " " deliminator
    keyword=split( readline( handle_in ) )

    # Check that the first column signal is "CRYST1"
    if keyword[1] == "CRYST1"
        # If so, loop over dimensions
        for i=1:3
            # Columns 2-4 are cell lengths
            cell.length[i] = parse( Float64, split( keyword )[ i + 1 ] )
            # Columns 5-7 are cell angles
            cell.angles[i] = parse( Float64, split( keyword )[ i + 4 ] )
        end
    # If not, send an error message and return false, false
    else
        print("Problem with .pdb file at: ",file_path,"\n")
        return false, false
    end

    # Initialize AtomMolList for atomic informations
    atoms = atom_mod.AtomMolList( nb_atoms )

    # Loop over atoms
    for atom=1:nb_atoms
        # Read line and parse with " " deliminator
        keyword = split( readline( handle_in ) )

        # Check that first column signal is "ATOM"
        if keyword[1] == "ATOM"
            # Parse 2nd column as atom index
            atoms.atom_index[atom] = parse( Int, keyword[2] )

            # Parse 3rd column as atom name
            atoms.atom_names[atom] = keyword[3]

            # Parse 4th column as molecule name
            atoms.mol_names[atom] =  keyword[4]

            # Parse 6th column as molecule index
            atoms.mol_index[atom] = parse( Int, keyword[5] )

            # Loop over dimensions
            for i=1:3
                # Parse columns 6-8 as atomic positions
                atoms.positions[ i, atom ] = parse( Float64, keyword[ 5 + i ] )
            end
        # If not, sends error message and return false, false
        else
            print("Problem with .pdb file at: ",file_path," at ATOM keyword.\n")
            return false, false
        end
    end

    # Close input file
    close( handle_in )

    #
    return atoms, cell
end
# Reads a .pdb file containing a trajectory with fixed cell
function readTrajAtomListFixedCell( file_path::T1 ) where { T1 <: AbstractString }
    # Argument
    # - file_path: path to the *.pdb file
    # Output
    # - traj: vector of AtomList, with trajectories of atoms
    # - cells: vector of Cell_param, with trajectory of cell

    # Get number of step in the trajectory
    nb_step = getNbSteps( file_path )

    # Get number of atoms in the traj
    nb_atoms = getNbAtoms( file_path )

    # Opens input file
    handle_in = open( file_path )

    # Initialize Cell_param for cell information
    cell = cell_mod.Cell_param()

    # Read and parse first line
    keyword = split( readline( handle_in ) )

    # Check that the first column is "CRYST1"
    if keyword[1] == "CRYST1"
        # Loop over dimensions
        for i=1:3
            # Columns 2-4 are cell lengths
            cell.length[ i ] = parse( Float64, keyword[ i + 1 ] )

            # Columns 5-7 are cell lengths
            cell.angles[ i ] = parse( Float64, keyword[ i + 4 ] )
        end
    # If not, sends a message and returns false, false
    else
        print("Problem with .pdb file at: ",file_path," at CRYST1 (start)\n")
        return false, false
    end

    # Goes back to the begining of the file
    seekstart( handle_in )

    # Initialize vector of AtomList for atomic trajectory
    traj = Vector{ AtomList }( undef, nb_step )

    # Loop over steps
    for step=1:nb_step
        # Read line and parse with " " deliminator
        keyword = split( readline( handle_in ) )

        # Check that the first column of line is "CRYST1"
        if keyword[1] != "CRYST1"
            # If not, returns false, false and sends error message
            print( "Problem with .pdb file at: ", file_path, " at CRYST1 step:", step, " !\n" )
            return false, false
        end

        # Check dummy
        check = false

        # Initialize ATomList for current step
        traj[step] = atom_mod.AtomList(nb_atoms)

        # Counter of atoms
        count_ = 1

        # Loop as long as we don't reach "END" keyword
        while ! check
            # Reads and parse line with " " deliminator
            keyword = split( readline( handle_in ) )

            # Check that keyword is "ATOM"
            if keyword[1] == "ATOM"
                # Get the atom name as third column
                traj[step].names = keyword[3]

                # Get atom counter as atom index
                traj[step].index = count_

                # Loop over dimension
                for i=1:3
                    # Columns 6-8 are parsed as atomic positions
                    traj[ step ].positions[ i, count_ ] = parse( Float64, keyword[ 5 + i ] )
                end

                # Increment atom counter
                count_ += 1

            # Else, check that the keyword is not "END", if so, stops steps
            elseif keyword[1] == "END"
                break
            end
        end
    end

    # Returns trajectory for atoms (traj) and cells (cell)
    return traj, cell
end
# Reads *.pdb file trajectory as vector of AtomList and Cell_param
function readTrajAtomList( file_path::T1 ) where { T1 <: AbstractString }
    # Argument
    # - file_path: path to the input file
    # Output
    # - traj: vector of AtomList with atomic positions
    # - cells: vector of Cell_param with cell informations

    # Get number of steps
    nb_step = getNbSteps( file_path )

    # Get number of atoms
    nb_atoms = getNbAtoms( file_path )

    # Open input file
    handle_in = open( file_path )

    # Initialize vector of AtomList for atom information
    traj = Vector{ AtomList }( undef, nb_step )

    # Initialize vector of Cell_param for cell information
    cells = Vector{cell_mod.Cell_param}(undef, nb_step )

    # Loop over steps
    for step=1:nb_step
        # Dummy to check for "END" signal
        check = false

        # Initialize local AtomList
        traj[step] = atom_mod.AtomList( nb_atoms )

        # Initialize atom counter
        count_ = 1

        # Loop as long as "END" signal is not reached
        while ! check
            # Read and parse line " " deliminator
            keyword = split( readline( handle_in ) )

            # Check that first column is CRYST1
            if keyword[1] == "CRYST1"
                # Initialize vectors for cell lengths and angles
                lengths = zeros(Real, 3 )
                angles  = zeros(Real, 3 )

                # Loop over dimensions
                for i=1:3
                    # Columns 2-4 are parsed as cell lenghts
                    lengths[i] = parse( Float64, keyword[ i + 1 ] )

                    # Columns 5-7 are parsed as cell lenghts
                    angles[i] = parse( Float64, keyword[ i + 4 ] )
                end

                # Compound information into Cell_param
                cells[ step ] = cell_mod.Cell_param( lengths, angles )
            # Check that the first column is "ATOM"
            elseif keyword[1] == "ATOM"
                # Parse second element as atom index
                traj[step].index[count_] = parse( Int, keyword[2] )

                # Parse third column as atom name
                traj[step].names[count_] = keyword[3]

                # Loop over dimensions
                for i=1:3
                    # Parse elements 6-8 as atomic positions
                    traj[step].positions[ i, count_ ] = parse( Float64, keyword[ 5 + i ] )
                end

                # Increment atom counter
                count_ += 1
            # If first column is END, end the step
            elseif keyword[1] == "END"
                break
            end
        end
    end

    # Returns vector of AtomList for atoms, vector of Cell_param for cell information
    return traj, cells
end
#-------------------------------------------------------------------------------

# Element of write first elements
#-------------------------------------------------------------------------------
# Creates and writes the CRYST1 line with cell information about the structure to a *.pdb file
# very parametric to fit all future evolution of the format
function writeCRYST1( handle_out::T1, lengths::Vector{T2}, angles::Vector{T3},
    round_length::T4=2,
    round_angles::T5=2,
    space_group::T6="P 1",
    z_number::T7=1,
    end_a::T8=15,
    end_b::T9=24,
    end_c::T10=33,
    end_alpha::T11=40,
    end_beta::T12=47,
    end_gamma::T13=54,
    start_sg::T14=56,
    end_z::T15=70 ) where { T1 <: IO, T2 <: Real, T3 <: Real, T4 <: Int, T5 <: Int, T6 <: AbstractString, T7 <: Int, T8 <: Int, T9 <: Int, T10 <: Int, T11 <: Int, T12 <: Int, T13 <: Int, T14 <: Int, T15 <: Int }
    # Arguments
    # - handle_out: IO stream to target file
    # - lengths:   lengths parameters of the cell
    # - angles:    angles parameters of the cell
    # Optionnal Arguments
    # - round_length: number of digits for lengths rounding
    # - round_angles: number of digits for angles rounding
    # - space_group: Space group in Hermann Mauguin notation, by default P 1 works...
    #              never seen a case where something else is used 28/01/2020 (M Moog)
    # - z_number: Number of polumeric chains or in hetero polymers
    #           the n umber of occurences of the most populous chain
    # - end_a:  Maximum colum for a (justify right for a )
    # - end_b:  Maximum colum for b (justify right for b )
    # - end_c:  Maximum colum for c (justify right for b )
    # - end_alpha:  Maximum colum for alpha (justify right for b )
    # - end_beta:  Maximum colum for beta (justify right for b )
    # - end_gamma: Maximum colum for gamma (justify right for b )
    # - start_sg: Min column for Space Group
    # - end_z: Maximum colum for gamma (justify right for b )
    # Output
    # - Bool: whether all went ok

    # Rounds cell lengths to the given digit
    a = string( round( lengths[1], digits=round_length ) )
    b = string( round( lengths[2], digits=round_length ) )
    c = string( round( lengths[3], digits=round_length ) )

    # Round cell angles to the given digits
    alpha = string( round( angles[1], digits=round_angles ) )
    beta  = string( round( angles[2], digits=round_angles ) )
    gamma = string( round( angles[3], digits=round_angles ) )

    # Casts Z number into a string
    z_number = string( z_number )

    # Writting CRYST1 signal
    cryst1 = string("CRYST1")

    # Adds information to the string with info about a into the appropriate columns number (right justified)
    cryst1 = addJustifyRight( end_a, cryst1, a )

    # If Something went wrong, sends error message and false
    if cryst1 == false
        print("Error writting CRYST1 line (a).\n")
        return false
    end

    # Adds information to the string with info about b into the appropriate columns number (right justified)
    cryst1 = addJustifyRight( end_b, cryst1, b )

    # If Something went wrong, sends error message and false
    if cryst1 == false
        print("Error writting CRYST1 line (b).\n")
        return false
    end

    # Adds information to the string with info about alpha into the appropriate columns number (right justified)
    cryst1 = addJustifyRight( end_c, cryst1, c )

    # If Something went wrong, sends error message and false
    if cryst1 == false
        print("Error writting CRYST1 line (c).\n")
        return false
    end

    # Adds information to the string with info about beta into the appropriate columns number (right justified)
    cryst1 = addJustifyRight( end_alpha, cryst1, alpha )

    # If Something went wrong, sends error message and false
    if cryst1 == false
        print("Error writting CRYST1 line (alpha).\n")
        return false
    end

    # Adds information to the string with info about gamma into the appropriate columns number (right justified)
    cryst1 = addJustifyRight( end_beta,  cryst1, beta  )

    # If Something went wrong, sends error message and false
    if cryst1 == false
        print("Error writting CRYST1 line (beta).\n")
        return false
    end

    # Adds information about the gamma factor to the string (right justified)
    cryst1 = addJustifyRight( end_gamma, cryst1, gamma )

    # If Something went wrong, sends error message and false
    if cryst1 == false
        print("Error writting CRYST1 line (gamma).\n")
        return false
    end

    # Adds information about the Space Group to the string (left justified)
    cryst1 = addJustifyLeft( start_sg, cryst1, space_group )

    # If Something went wrong, sends error message and false
    if cryst1 == false
        print("Error writting CRYST1 line (Space Group).\n")
        return false
    end

    # Z number (right justified)
    cryst1 = addJustifyRight( end_z, cryst1, z_number )

    # If Something went wrong, sends error message and false
    if cryst1 == false
        print("Error writting CRYST1 line (z).\n")
        return false
    end

    # Adds line break at the end of the line
    cryst1 = string( cryst1, "\n" )

    # Writting line
    Base.write( handle_out, cryst1 )

    # Return true if all went ok
    return true
end
# Creates and writes the MODEL line with a name of structure to *.pdb file
# NB: not the most useful functions ever but if it works...
function writeMODEL( handle_out::T1, name_model::T2 )  where { T1 <: IO, T2 <: AbstractString }
    # Argument
    # - handle_out: IO handler for the output file
    # - name_model: name of the model to write
    # Output
    # - Bool: whether it worked or not

    # Writting to the file
    write( handle_out, string("MODEL ",name_model))

    # Returning true
    return true
end
# Creates and writes the ATOM line with atomic information about a single atom
# Parametric, works if there is or isn't molecules
function writeATOM( handle_out::T1, atom_index::T2, atom_name::T3, atom_positions::Vector{T4},
    max_col_atom_index::T5=11,
    max_col_atom_name::T26=14,
    max_col_alt_loc::T6=17,
    max_col_res_name::T7=20,
    max_col_mol_name::T8=22,
    max_col_mol_index::T9=26,
    max_col_ins_res_code::T10=27,
    max_col_x::T11=38,
    max_col_y::T12=46,
    max_col_z::T13=54,
    max_col_occup::T14=60,
    max_col_temp_fac::T15=66,
    max_col_atom_name2::T16=78,
    max_col_atom_charge::T17=80,
    alt_location::T18=" ",
    residue_name::T19="   ",
    molecule_name::T20="X",
    molecule_index::T21=0,
    residue_insertion_code::T22=" ",
    occupancy::T23=0,
    tempfac::T24=0,
    charge::T25=0 ) where { T1 <: IO, T2 <: Int, T3 <: AbstractString, T4 <: Int, T5 <: Int, T6 <: Int, T7 <: Int,
    T8 <: Int, T9 <: Int, T10 <: Int, T11 <: Int, T12 <: Int, T13 <: Int, T14 <: Int, T15 <: Int, T16 <: Int, T17 <: Int,
    T18 <: AbstractString, T19 <: AbstractString, T20 <: AbstractString, T21 <: Int, T22 <: AbstractString, T23 <: Real,
    T24 <: Real, T25 <: Real, T26 <: Int }
    # Argument
    # - handle_out: IO for the output file
    # - atom_name: name of the target atom to add to the file
    # - atom_index: index of the target atom to add to the file
    # - atom_positions: positions of the target atom to add to the file
    # Optionnal Argument
    # - max_col_atom_index: maximum column for atom index
    # - max_col_alt_loc: maximum column for alternative location
    # - max_col_res_name: maximum column for residue number
    # - max_col_mol_name: maximum column for molecule name
    # - max_col_mol_index: maximum column for molecule index
    # - max_col_ins_res_code: maximum column for insertion code of residue
    # - max_col_x: maximum column for x position
    # - max_col_y: maximum column for y position
    # - max_col_z: maximum column for z position
    # - max_col_occup: maximum column for occupancy of atom
    # - max_col_temp_fac: maximum column for temperature factor
    # - max_col_atom_name2: maximum column for second atom name
    # - max_col_atom_charge: maximum column for atomic charge
    # - alt_location: alternative location signal
    # - residue_name: name of the residue to which the atom belong
    # - molecule_name: name of the molecule to which the atom belong
    # - molecule_index: index of the molecule to which the atom belong
    # - residue_insertion_code: code of atom insertion
    # - occupancy: occupancy of the atom
    # - tempfac: temperature factor
    # - charge: charge of the atom
    # Output
    # - Bool: whether it went ok

    # Start with "ATOM" keyword
    atom_line = string("ATOM")

    # Adds the index of atom
    atom_line = utils.addJustifyRight( max_col_atom_index, string(atom_index), atom_line )

    # Adds the name of atom
    atom_line = utils.addJustifyRight( max_col_atom_name, atom_name, atom_line )

    # Alternative location symbol
    atom_line = utils.addJustifyRight( max_col_alt_loc, alt_location, atom_line )

    # Adds the name of the residue (for protein)
    atom_line = utils.addJustifyRight( max_col_res_name, residue_name, atom_line )

    # Adds the name of the chain/molecule
    atom_line = utils.addJustifyRight( max_col_mol_name, molecule_name, atom_line )

    # Adds the index of the molecule
    atom_line = utils.addJustifyRight( max_col_mol_index, string( molecule_index ), atom_line )

    # Adds the insertion code of the residue
    atom_line = utils.addJustifyRight( max_col_ins_res_code, residue_insertion_code, atom_line )

    # Adds the position x
    atom_line = utils.addJustifyRight( max_col_x, string( atom_positions[1] ), atom_line )

    # Adds the position y
    atom_line = utils.addJustifyRight( max_col_y, string( atom_positions[2] ), atom_line )

    # Adds the position z
    atom_line = utils.addJustifyRight( max_col_z, string( atom_positions[3] ), atom_line )

    # Adds the occupancy
    atom_line = utils.addJustifyRight( max_col_occup, string( occupancy ), atom_line )

    # Adds the temperature factor
    atom_line = utils.addJustifyRight( max_col_temp_fac, string( tempfac ), atom_line )

    # Adds the name of the atom, again
    atom_line = utils.addJustifyRight( max_col_atom_name2, atom_name, atom_line )

    # Adds the charge
    atom_line = utils.addJustifyRight( max_col_atom_charge, string( charge ), atom_line )

    # Adds the end line to the string
    atom_line = string( atom_line, "\n")

    # If something went wrong, returns false
    if atom_line == false
        return false
    end

    # If we reach this point, we return true
    return true
end
#-------------------------------------------------------------------------------

# Writing file
#-------------------------------------------------------------------------------
# Writes a PDB file according to standard format (2011) for a single atomic structure in AtomList, with a Cell_param for cell
function writePdb( file_path::T1, atoms::T2, cell::T3 ) where { T1 <: AbstractString , T2 <: atom_mod.AtomList, T3 <: cell_mod.Cell_param }
    # Argument
    # - file_path: path to the input file
    # - atoms: AtomList with atomic information
    # - cell: Cell_param with cell informations
    # Output
    # - Bool: whether it worked or not

    # Opens output file
    handle_out = open( file_path, "w" )

    # Writes CRYST1 line with cell information
    if ! writeCRYST1( handle_out, cell.length, cell.angles )
        # Returns false if something went wrong
        return false
    end

    # Get number of atoms in the structure
    nb_atoms = size(atoms.names)[1]

    # Loop over atoms
    for atom=1:nb_atoms
        writeATOM( handle_out, atoms.index[ atom ], atoms.name[ atom ], atoms.positions[ :, atom ] )
    end

    # Writes "END" signal for the step
    Base.write( handle_out, "END\n" )

    # Close output file
    close( handle_out )

    # Returns true if something went ok
    return true
end
# Writes a PDB file according to standard format (2011) for several atomic structure (vector of AtomList) and a single cell (Cell_param)
function writePdb( file_path::T1, traj::Vector{T2}, cell::T3 ) where { T1 <: AbstractString , T2 <: atom_mod.AtomList, T3 <: cell_mod.Cell_param }
    # Argument
    # - file_path: path to the output file
    # - traj: vector of AtomList with atomic trajectory
    # - cell: Cell_param with cell information
    # Output
    # - None

    # Open output file
    handle_out = open( file_path, "w" )

    # Get number of step of the trajectory
    nb_step = size(traj)

    # Loop over steps
    for step=1:nb_step

        if ! writeCRYST1( handle_out, cell.length, cell.angles )
            print("Error writting cell information!\n")
            return false
        end

        # Get number of atoms
        nb_atoms = size(traj[step].names)[1]

        for atom=1:nb_atoms
            writeATOM( handle_out, traj[ step ].index[ atom ], traj[ step ].name[ atom ], traj[ step ].positions[ :, atom ] )
        end

        Base.write(handle_out,"END\n")
    end

    close(handle_out)

    return
end
# Writes a PDB file according to standard format (2011) for several atomic structure (vector of AtomList) and cells (vector of Cell_param)
function writePdb( file_path::T1, traj::Vector{T2}, cells::Vector{T3} ) where { T1 <: AbstractString , T2 <: atom_mod.AtomList, T3 <: cell_mod.Cell_param }
    # Argument
    # - file_path: path to the output file
    # - traj: vector of Atomlist with atomic information for the trajectory
    # - cells: vector of Cell_param with cell information for the trajectory
    # Output
    # - Bool; whether writting was successful

    # Open output file
    handle_out = open( file_path, "w" )

    # Get number of steps
    nb_step=size(traj)[1]

    # Loop over steps
    for step=1:nb_step

        # Write cell information
        if ! writeCRYST1( handle_out, cells[step].length, cells[step].angles )
            print("Error writting cell information!\n")
            return false
        end

        # Atomic Positions
        nb_atoms = size(traj[step].names)[1]

        # Loop over atoms
        for atom=1:nb_atoms
            # Write line for current atom
            if ! writeATOM( handle_out, traj[ step ].index[ atom ], traj[ step ].name[ atom ], traj[ step ].positions[ :, atom ] )
                print("Error writting atom: ", atom, ".\n")
                return false
            end
        end

        # Write END signal for step
        Base.write(handle_out,"END\n")
    end

    # Close output file
    close(handle_out)

    # If we reach this point returns true as all went well
    return true
end
# Writes a PDB file according to standard format (2011) for several atomic structure (vector of AtomList) and cells (vector of Cell_param) using IO handle for file
function writePdb( handle_out::T1, traj::Vector{T2}, cells::Vector{T3} ) where { T1 <: IO , T2 <: atom_mod.AtomList, T3 <: cell_mod.Cell_param }
    # Argument
    # - handle_out: IO handler for output file
    # - traj: vector of AtomList, contains atoms trajectory
    # - cells: vector of Cell_param, contains cell trajectory
    # Output
    # Bool: whether writting went ok

    # Get number of steps
    nb_step = size( traj )[1]

    # Loop over steps
    for step=1:nb_step
        # Write cell information
        if ! writeCRYST1( handle_out, cells[step].length, cells[step].angles )
            print("Error writting cell information!\n")
            return false
        end

        # Atomic Positions
        nb_atoms = size(traj[step].names)[1]

        # Loop over atoms
        for atom=1:nb_atoms
            # Write atoms to file
            writeATOM( handle_out, traj[ step ].index[ atom ], traj[ step ].name[ atom ], traj[ step ].positions[ :, atom ] )
        end

        # Write END of line signal
        Base.write(handle_out,"END\n")
    end

    # Returns true if we reach this point
    return true
end
#
function writePdbPivClustering( file::T1, atoms::T2, cell::T3 ) where { T1 <: AbstractString , T2 <: atom_mod.AtomList, T3 <: cell_mod.Cell_param }

    # Writes a PDB file according to standard format (2011)

    handle_out=open(file,"w")

    if ! writeCRYST1( handle_out, cell.length, cell.angles )
        print("Error writting cell information!\n")
        return false
    end

    if ! writeMODEL( handle_out, "X")
        print("Error writting MODEL line !\n")
        return false
    end

    # Atomic Positions
    nb_atoms = size(atoms.names)[1]
    for i=1:nb_atoms
        atom="ATOM"
        # Right justified
        atom=utils.spaces(atom,11-length(string(atoms.index[i]))-length(atom))
        atom=string(atom,atoms.index[i])
        # If atom name is 1 length, start on 13, otherwise 14
        # atom names are left justified here
        if length(atoms.names[i])==1
            atom=utils.spaces(atom,13-length(atom))
        else
            atom=utils.spaces(atom,14-length(atom))
        end
        atom=string(atom,atoms.names[i])

        # Alternate location indicator ( default blank )
        alt_loc=string(" ")
        atom=utils.spaces(atom,17-length(atom)-length(alt_loc))
        atom=string(atom," ")

        # Residue name (3 col max) righ justified?
        residue_name=string("   ")
        atom=utils.spaces(atom,20-length(atom)-length(residue_name))
        atom=string(atom,residue_name)

        # Chain Identifier
        # Default is "X"
        chain_id=string("X")
        atom=utils.spaces(atom,22-length(atom)-length(chain_id))
        atom=string(atom,"X")

        # Residue Sequence Nb - Molecule nb
        mol_nb=string("1")
        atom=utils.spaces(atom,26-length(atom)-length(mol_nb))
        atom=string(atom,mol_nb)

        # Code for insertion of residues
        code_insertion=string(" ")
        atom=utils.spaces(atom,27-length(atom)-length(code_insertion))
        atom=string(atom,code_insertion)

        # Positions are right justified
        # X
        x=string( round(atoms.positions[i,1], digits=3 ) )
        atom=utils.spaces( atom, 38-length(atom)-length(x) )
        atom=string( atom, x)
        # Y
        y=string( round(atoms.positions[i,2], digits=3 ) )
        atom=utils.spaces( atom, 46-length(atom)-length(y) )
        atom=string( atom, y )
        # Z
        z=string( round( atoms.positions[i,3], digits=3 ))
        atom=utils.spaces( atom, 54-length(atom)-length(z) )
        atom=string( atom, z )

        # Occupancy and Temperature Factor (default 0, except for PLUMED where used for masses)
        # Occupancy
        occ = periodicTable.names2Z( atoms.names[i] )
        atom=utils.spaces(atom, 60-length(atom)-length(occ) )
        atom=string(atom, occ )
        # Temperature Factor
        tempfac = periodicTable.names2Z( atoms.names[i] )
        atom=utils.spaces(atom,66-length(atom))
        atom=string(atom, string( 0.0 ),string(0) )

        # Atom name (right justified)
        atom=utils.spaces(atom,78-length(atom)-length(atoms.names[i]))
        atom=string(atom,atoms.names[i])

        # Charge (2 col max, default is blank)
        charge=string(" ")
        atom=utils.spaces(atom,80-length(atom)-length(charge))
        atom=string(atom,charge)

        # End of line
        atom=string(atom,"\n")
        Base.write(handle_out,atom)
    end

    Base.write(handle_out,"END\n")

    close(handle_out)

    return
end
#
function writePdbPivClustering( file::T1, traj::Vector{T2}, cells::Vector{T3} ) where { T1 <: AbstractString , T2 <: atom_mod.AtomList, T3 <: cell_mod.Cell_param }

    nb_step  = size( traj )[1]
    nb_atoms = size( traj[1].names )[1]

    # Writes a PDB file according to standard format (2011)
    handle_out=open(file,"w")

    for step=1:nb_step

        if ! writeCRYST1( handle_out, cells[step].length, cells[step].angles )
            print("Error writting cell information!\n")
            return false
        end

        if ! writeMODEL( handle_out, "X")
            print("Error writting MODEL line !\n")
            return false
        end

        # Atomic Positions
        for i=1:nb_atoms
            atom="ATOM"
            # Right justified
            atom = utils.spaces( atom, 11 - length( string(traj[step].index[i] ) ) - length( atom ) )
            atom = string( atom, traj[step].index[i] )
            # If atom name is 1 length, start on 13, otherwise 14
            # atom names are left justified here
            if length( traj[step].names[i] ) == 1
                atom = utils.spaces( atom, 13 - length( atom ) )
            else
                atom = utils.spaces( atom, 14 - length( atom ) )
            end
            atom = string( atom, traj[step].names[i] )

            # Alternate location indicator ( default blank )
            alt_loc = string( " " )
            atom = utils.spaces( atom, 17 - length(atom) - length( alt_loc ) )
            atom = string( atom, " " )

            # Residue name (3 col max) righ justified?
            residue_name = string( "   " )
            atom = utils.spaces( atom, 20 - length(atom) - length( residue_name ) )
            atom = string( atom, residue_name )

            # Chain Identifier
            # Default is "X"
            chain_id = string( "X" )
            atom     = utils.spaces( atom, 22 - length( atom ) - length( chain_id ) )
            atom     = string( atom, "X" )

            # Residue Sequence Nb - Molecule nb
            mol_nb = string( "1" )
            atom = utils.spaces( atom, 26 - length( atom ) - length( mol_nb ) )
            atom = string( atom, mol_nb )

            # Code for insertion of residues
            code_insertion=string(" ")
            atom=utils.spaces(atom,27-length(atom)-length(code_insertion))
            atom=string(atom,code_insertion)

            # Positions are right justified
            # X
            x    = string( round( traj[step].positions[i,1], digits=3 ) )
            atom = utils.spaces( atom, 38 - length(atom) - length(x) )
            atom = string( atom, x)
            # Y
            y    = string( round( traj[step].positions[i,2], digits=3 ) )
            atom = utils.spaces( atom, 46 - length(atom) - length(y) )
            atom = string( atom, y )
            # Z
            z    = string( round( traj[step].positions[i,3], digits=3 ) )
            atom = utils.spaces( atom, 54 - length(atom) - length(z) )
            atom = string( atom, z )

            # Occupancy and Temperature Factor (default 0, except for PLUMED where used for masses)
            # Occupancy
            occ =  periodicTable.names2Z( traj[step].names[i] )
            atom = utils.spaces( atom, 60 - length(atom) - length(occ) )
            atom = string(atom, occ )
            # Temperature Factor
            tempfac = periodicTable.names2Z( traj[step].names[i] )
            atom    = utils.spaces( atom, 66 - length(atom) )
            atom    = string( atom, string( 0.0 ), string(0) )

            # Atom name (right justified)
            atom = utils.spaces( atom, 78 - length(atom) - length(traj[step].names[i]) )
            atom = string( atom, traj[step].names[i] )

            # Charge (2 col max, default is blank)
            charge = string( " " )
            atom   = utils.spaces( atom, 80 - length(atom) - length(charge) )
            atom   = string( atom, charge )

            # End of line
            atom = string( atom, "\n" )
            Base.write( handle_out, atom )
        end

        Base.write( handle_out, "END\n" )

    end

    close(handle_out)

    return true
end
#
function writePdbPlumed( atoms::T1, cell::T2, file::T3 ) where { T1 <: atom_mod.AtomMolList, T2 <: cell_mod.Cell_param, T3 <: AbstractString }

  out=open(file,"w")

  a,b,c = string(cell.length[1]), string(cell.length[2]), string(cell.length[3])
  alpha, beta, gamma = string(cell.angles[1]), string(cell.angles[2]), string(cell.angles[3])

  cryst1=string("CRYST1 ",a)
  cryst1=utils.spaces(cryst1,16-length(cryst1))
  cryst1=string(cryst1,b)
  cryst1=utils.spaces(cryst1,25-length(cryst1))
  cryst1=string(cryst1,c)
  cryst1=utils.spaces(cryst1,34-length(cryst1))
  cryst1=string(cryst1,alpha)
  cryst1=utils.spaces(cryst1,41-length(cryst1))
  cryst1=string(cryst1,beta)
  cryst1=utils.spaces(cryst1,48-length(cryst1))
  cryst1=string(cryst1,gamma)
  cryst1=utils.spaces(cryst1,56-length(cryst1))
  cryst1=string(cryst1,"P 1")
  cryst1=utils.spaces(cryst1,67-length(cryst1))
  cryst1=string(cryst1,"1")
  cryst1=string(cryst1,"\n")
  Base.write(out,cryst1)

  nb_atoms = size(atoms.atom_names)[1]
  for i=1:nb_atoms
    atom="ATOM"
    atom=utils.spaces(atom,11-length(string(atoms.atom_index[i]))-length(atom))
    atom=string(atom,atoms.atom_index[i])
    atom=utils.spaces(atom,13-length(atom))
    atom=string(atom,atoms.atom_names[i])
    atom=utils.spaces(atom,17-length(atom))
    atom=string(atom,atoms.mol_names[i])
    atom=utils.spaces(atom,21-length(atom))
    atom=string(atom,"X")
    atom=utils.spaces(atom,25-length(atom))
    atom=string(atom,atoms.mol_index[i])
    atom=utils.spaces(atom,32-length(atom))
    atom=string(atom,round(atoms.positions[i,1],digits=3 ) )
    atom=utils.spaces(atom,40-length(atom))
    atom=string(atom,round(atoms.positions[i,2],digits=3 ) )
    atom=utils.spaces(atom,48-length(atom))
    atom=string(atom,round(atoms.positions[i,3],digits=3 ) )
    atom=utils.spaces(atom,56-length(atom))
    atom=string(atom, string( round( periodicTable.names2Z(atoms.names[i]) ) ) )
    atom=utils.spaces(atom,62-length(atom))
    atom=string(atom, string( round( periodicTable.names2Z(atoms.names[i]) ) ) )
    atom=utils.spaces(atom,77-length(atom))
    atom=string(atom,"\n")
    Base.write(out,atom)
  end

  Base.write(out,"END\n")

  close(out)

  return
end
#
function writePdbPlumed( atoms::T1, cell::T2, file::T3 ) where { T1 <: atom_mod.AtomList, T2 <: cell_mod.Cell_param, T3 <: AbstractString }

  out=open(file,"w")

  a,b,c = string(round(cell.length[1],digits=2)), string(round(cell.length[2],digits=2)), string(round(cell.length[3],digits=2))
  alpha, beta, gamma = string(round(cell.angles[1],digits=2)), string(round(cell.angles[2],digits=2)), string(round(cell.angles[3],digits=2))

  cryst1=string("CRYST1 ",a)
  cryst1=utils.spaces(cryst1,16-length(cryst1))
  cryst1=string(cryst1,b)
  cryst1=utils.spaces(cryst1,25-length(cryst1))
  cryst1=string(cryst1,c)
  cryst1=utils.spaces(cryst1,34-length(cryst1))
  cryst1=string(cryst1,alpha)
  cryst1=utils.spaces(cryst1,41-length(cryst1))
  cryst1=string(cryst1,beta )
  cryst1=utils.spaces(cryst1,48-length(cryst1))
  cryst1=string(cryst1,gamma )
  cryst1=utils.spaces(cryst1,56-length(cryst1))
  cryst1=string(cryst1,"P 1")
  cryst1=utils.spaces(cryst1,67-length(cryst1))
  cryst1=string(cryst1,"1")
  cryst1=string(cryst1,"\n")
  Base.write(out,cryst1)

  nb_atoms = size(atoms.names)[1]
  for i=1:nb_atoms
    atom="ATOM"
    atom=utils.spaces(atom,11-length(string(atoms.index[i]))-length(atom))
    atom=string(atom,atoms.index[i])
    atom=utils.spaces(atom,13-length(atom))
    atom=string(atom,atoms.names[i])
    atom=utils.spaces(atom,17-length(atom))
    atom=string(atom," ")
    atom=utils.spaces(atom,21-length(atom))
    atom=string(atom,"X")
    atom=utils.spaces(atom,25-length(atom))
    atom=string(atom,1)
    atom=utils.spaces(atom,32-length(atom))
    atom=string(atom,round(atoms.positions[i,1], digits=3 ) )
    atom=utils.spaces(atom,40-length(atom))
    atom=string(atom,round(atoms.positions[i,2], digits=3 ) )
    atom=utils.spaces(atom,48-length(atom))
    atom=string(atom,round(atoms.positions[i,3], digits=3 ) )
    atom=utils.spaces(atom,56-length(atom))
    atom=string(atom, string( 0.0 ),string(0) )
    atom=utils.spaces(atom,62-length(atom))
    atom=string(atom, string( 0.0 ),string(0) )
    atom=utils.spaces(atom,76-length(atom))
    atom=string(atom,atoms.names[i])
    atom=string(atom,"\n")
    Base.write(out,atom)
  end

  Base.write(out,"END\n")

  close(out)

  return
end
#-------------------------------------------------------------------------------
end
