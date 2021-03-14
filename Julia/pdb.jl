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

# Writting file in justification file
#-------------------------------------------------------------------------------
function addJustifyRight( max_column::T1, string_target::T2, string_toadd::T3 ) where  { T1 <: Int, T2 <: AbstractString, T3 <: AbstractString }
    string_target = utils.spaces( string_target, max_column-length(string_target)-length(string_toadd) )
    return string( string_target, string_toadd)
end
function addJustifyLeft( max_column::T1, string_target::T2, string_toadd::T3 ) where { T1 <: Int, T2 <: AbstractString, T3 <: AbstractString }
    string_target = utils.spaces( string_target, max_column-length(string_target) )
    return string( string_target, string_toadd)
end
#-------------------------------------------------------------------------------

# Element of write first elements
#-------------------------------------------------------------------------------
function writeCRYST1(  handle_out::IO, lengths::Vector{Real}, angles::Vector{Real},
    round_length::Int=2,
    round_angles::Int=2,
    space_group::AbstractString="P 1",
    z_number::Int=1,
    end_a::Int=15,
    end_b::Int=24,
    end_c::Int=33,
    end_alpha::Int=40,
    end_beta::Int=47,
    end_gamma::Int=54,
    start_sg::Int=56,
    end_z::Int=70 )

    # Role: write CRYST1 to file
    # Parametric, in order to adapt to the various PDB format in existence

    # Arguments
    # handle_out: IO stream to target file
    # lengths:   lengths parameters of the cell
    # angles:    angles parameters of the cell

    # Optionnal Arguments
    # round_length: number of digits for lengths rounding
    # round_angles: number of digits for angles rounding
    # space_group: Space group in Hermann Mauguin notation, by default P 1 works...
    #              never seen a case where something else is used 28/01/2020 (M Moog)
    # z_number: Number of polumeric chains or in hetero polymers
    #           the n umber of occurences of the most populous chain
    # end_a:  Maximum colum for a (justify right for a )
    # end_b;  Maximum colum for b (justify right for b )
    # end_c;  Maximum colum for c (justify right for b )
    # end_alpha;  Maximum colum for alpha (justify right for b )
    # end_beta;  Maximum colum for beta (justify right for b )
    # end_gamma;  Maximum colum for gamma (justify right for b )
    # start_sg;  Min column for Space Group
    # end_z;  Maximum colum for gamma (justify right for b )

    # lengths
    a = string( round( lengths[1], digits=round_length ) )
    b = string( round( lengths[2], digits=round_length ) )
    c = string( round( lengths[3], digits=round_length ) )
    # angles
    alpha = string( round( angles[1], digits=round_angles ) )
    beta  = string( round( angles[2], digits=round_angles ) )
    gamma = string( round( angles[3], digits=round_angles ) )
    # z number
    z_number=string(z_number)

    # Writting CRYST1 info -> Cell information
    cryst1 = string("CRYST1")
    # Cell lengths (right justified)
    cryst1 = addJustifyRight( end_a, cryst1, a )
    if cryst1 == false
        print("Error writting CRYST1 line (a).\n")
        return false
    end
    cryst1 = addJustifyRight( end_b, cryst1, b )
    if cryst1 == false
        print("Error writting CRYST1 line (b).\n")
        return false
    end
    cryst1 = addJustifyRight( end_c, cryst1, c )
    if cryst1 == false
        print("Error writting CRYST1 line (c).\n")
        return false
    end
    # Angles (right justified)
    cryst1 = addJustifyRight( end_alpha, cryst1, alpha )
    if cryst1 == false
        print("Error writting CRYST1 line (alpha).\n")
        return false
    end
    cryst1 = addJustifyRight( end_beta,  cryst1, beta  )
    if cryst1 == false
        print("Error writting CRYST1 line (beta).\n")
        return false
    end
    cryst1 = addJustifyRight( end_gamma, cryst1, gamma )
    if cryst1 == false
        print("Error writting CRYST1 line (gamma).\n")
        return false
    end
    # Space Group (left justified)
    cryst1 = addJustifyLeft( start_sg, cryst1, space_group )
    if cryst1 == false
        print("Error writting CRYST1 line (Space Group).\n")
        return false
    end
    # Z number (right justified)
    cryst1 = addJustifyRight( end_z, cryst1, z_number )
    if cryst1 == false
        print("Error writting CRYST1 line (z).\n")
        return false
    end
    cryst1=string(cryst1,"\n")

    # Writting line
    Base.write(handle_out,cryst1)

    return true
end
function writeMODEL( handle_out::T1, name_model::T2 )  where { T1 <: IO, T2 <: AbstractString }
    write( handle_out, string("MODEL ",name_model))
    return true
end
#-------------------------------------------------------------------------------

# Writing file
#-------------------------------------------------------------------------------
function writePdb( file::T1, atoms::T2, cell::T3 ) where { T1 <: AbstractString , T2 <: atom_mod.AtomList, T3 <: cell_mod.Cell_param }

    # Writes a PDB file according to standard format (2011)

    handle_out=open(file,"w")

    if ! writeCRYST1( handle_out, cell.length, cell.angles )
        print("Error writting cell information!\n")
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
        if length(atoms.names[i]) == 1
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
function writePdb( file::T1, traj::Vector{T2}, cell::T3 ) where { T1 <: AbstractString , T2 <: atom_mod.AtomList, T3 <: cell_mod.Cell_param }

    # Writes a PDB file according to standard format (2011)

    handle_out=open(file,"w")

    nb_step=size(traj)
    for step=1:nb_step

        if ! writeCRYST1( handle_out, cell.length, cell.angles )
            print("Error writting cell information!\n")
            return false
        end

        # Atomic Positions
        nb_atoms = size(traj[step].names)[1]
        for i=1:nb_atoms
            atom="ATOM"
            # Right justified
            atom=utils.spaces(atom,11-length(string(traj[step].index[i]))-length(atom))
            atom=string(atom,traj[step].index[i])
            # If atom name is 1 length, start on 13, otherwise 14
            # atom names are left justified here
            if length(traj[step].names[i])==1
                atom=utils.spaces(atom,13-length(atom))
            else
                atom=utils.spaces(atom,14-length(atom))
            end
            atom=string(atom,traj[step].names[i])

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
            x=string( round(traj[step].positions[i,1], digits=3 ) )
            atom=utils.spaces( atom, 38-length(atom)-length(x) )
            atom=string( atom, x)
            # Y
            y=string( round(traj[step].positions[i,2], digits=3 ) )
            atom=utils.spaces( atom, 46-length(atom)-length(y) )
            atom=string( atom, y )
            # Z
            z=string( round( traj[step].positions[i,3], digits=3 ))
            atom=utils.spaces( atom, 54-length(atom)-length(z) )
            atom=string( atom, z )

            # Occupancy and Temperature Factor (default 0, except for PLUMED where used for masses)
            # Occupancy
            occ = periodicTable.names2Z( traj[step].names[i] )
            atom=utils.spaces(atom, 60-length(atom)-length(occ) )
            atom=string(atom, occ )
            # Temperature Factor
            tempfac = periodicTable.names2Z( traj[step].names[i] )
            atom=utils.spaces(atom,66-length(atom))
            atom=string(atom, string( 0.0 ),string(0) )

            # Atom name (right justified)
            atom=utils.spaces(atom,78-length(atom)-length(traj[step].names[i]))
            atom=string(atom,traj[step].names[i])

            # Charge (2 col max, default is blank)
            charge=string(" ")
            atom=utils.spaces(atom,80-length(atom)-length(charge))
            atom=string(atom,charge)

            # End of line
            atom=string(atom,"\n")
            Base.write(handle_out,atom)
        end

        Base.write(handle_out,"END\n")
    end

    close(handle_out)

    return
end
function writePdb( file::T1, traj::Vector{T2}, cells::Vector{T3} ) where { T1 <: AbstractString , T2 <: atom_mod.AtomList, T3 <: cell_mod.Cell_param }

    # Writes a PDB file according to standard format (2011)

    handle_out=open(file,"w")

    nb_step=size(traj)[1]
    for step=1:nb_step

        if ! writeCRYST1( handle_out, cells[step].length, cells[step].angles )
            print("Error writting cell information!\n")
            return false
        end

        # Atomic Positions
        nb_atoms = size(traj[step].names)[1]
        for i=1:nb_atoms
            atom="ATOM"
            # Right justified
            atom=utils.spaces(atom,11-length(string(traj[step].index[i]))-length(atom))
            atom=string(atom,traj[step].index[i])
            # If atom name is 1 length, start on 13, otherwise 14
            # atom names are left justified here
            if length(traj[step].names[i])==1
                atom=utils.spaces(atom,13-length(atom))
            else
                atom=utils.spaces(atom,14-length(atom))
            end
            atom=string(atom,traj[step].names[i])

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
            x=string( round(traj[step].positions[i,1], digits=3 ) )
            atom=utils.spaces( atom, 38-length(atom)-length(x) )
            atom=string( atom, x)
            # Y
            y=string( round(traj[step].positions[i,2], digits=3 ) )
            atom=utils.spaces( atom, 46-length(atom)-length(y) )
            atom=string( atom, y )
            # Z
            z=string( round( traj[step].positions[i,3], digits=3 ))
            atom=utils.spaces( atom, 54-length(atom)-length(z) )
            atom=string( atom, z )

            # Occupancy and Temperature Factor (default 0, except for PLUMED where used for masses)
            # Occupancy
            occ = periodicTable.names2Z( traj[step].names[i] )
            atom=utils.spaces(atom, 60-length(atom)-length(occ) )
            atom=string(atom, occ )
            # Temperature Factor
            tempfac = periodicTable.names2Z( traj[step].names[i] )
            atom=utils.spaces(atom,66-length(atom))
            atom=string(atom, string( 0.0 ),string(0) )

            # Atom name (right justified)
            atom=utils.spaces(atom,78-length(atom)-length(traj[step].names[i]))
            atom=string(atom,traj[step].names[i])

            # Charge (2 col max, default is blank)
            charge=string(" ")
            atom=utils.spaces(atom,80-length(atom)-length(charge))
            atom=string(atom,charge)

            # End of line
            atom=string(atom,"\n")
            Base.write(handle_out,atom)
        end

        Base.write(handle_out,"END\n")
    end

    close(handle_out)

    return
end
function writePdb( handle_out::T1, traj::Vector{T2}, cells::Vector{T3} ) where { T1 <: IO , T2 <: atom_mod.AtomList, T3 <: cell_mod.Cell_param }

    # Writes a PDB file according to standard format (2011)

    nb_step=size(traj)[1]
    for step=1:nb_step

        if ! writeCRYST1( handle_out, cells[step].length, cells[step].angles )
            print("Error writting cell information!\n")
            return false
        end

        # Atomic Positions
        nb_atoms = size(traj[step].names)[1]
        for i=1:nb_atoms
            atom="ATOM"
            # Right justified
            atom=utils.spaces(atom,11-length(string(traj[step].index[i]))-length(atom))
            atom=string(atom,traj[step].index[i])
            # If atom name is 1 length, start on 13, otherwise 14
            # atom names are left justified here
            if length(traj[step].names[i])==1
                atom=utils.spaces(atom,13-length(atom))
            else
                atom=utils.spaces(atom,14-length(atom))
            end
            atom=string(atom,traj[step].names[i])

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
            x=string( round(traj[step].positions[i,1], digits=3 ) )
            atom=utils.spaces( atom, 38-length(atom)-length(x) )
            atom=string( atom, x)
            # Y
            y=string( round(traj[step].positions[i,2], digits=3 ) )
            atom=utils.spaces( atom, 46-length(atom)-length(y) )
            atom=string( atom, y )
            # Z
            z=string( round( traj[step].positions[i,3], digits=3 ))
            atom=utils.spaces( atom, 54-length(atom)-length(z) )
            atom=string( atom, z )

            # Occupancy and Temperature Factor (default 0, except for PLUMED where used for masses)
            # Occupancy
            occ = periodicTable.names2Z( traj[step].names[i] )
            atom=utils.spaces(atom, 60-length(atom)-length(occ) )
            atom=string(atom, occ )
            # Temperature Factor
            tempfac = periodicTable.names2Z( traj[step].names[i] )
            atom=utils.spaces(atom,66-length(atom))
            atom=string(atom, string( 0.0 ),string(0) )

            # Atom name (right justified)
            atom=utils.spaces(atom,78-length(atom)-length(traj[step].names[i]))
            atom=string(atom,traj[step].names[i])

            # Charge (2 col max, default is blank)
            charge=string(" ")
            atom=utils.spaces(atom,80-length(atom)-length(charge))
            atom=string(atom,charge)

            # End of line
            atom=string(atom,"\n")
            Base.write(handle_out,atom)
        end

        Base.write(handle_out,"END\n")
    end

    return true
end
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
