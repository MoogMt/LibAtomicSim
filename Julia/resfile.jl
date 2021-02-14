module resfile

# Loading necessary modules from LibAtomicSim
using utils
using conversion
using atom_mod
using cell_mod
using filexyz
using pdb

# Loading necessary module from general repository
using LinearAlgebra

# Cell information
#-------------------------------------------------------------------------------
# Read cell info from .res file using IO handler for the input file
function extractCellInfo( file_io::T1 ) where { T1 <: IO}
    # Argument
    # - file_io: handler of the input .res file
    # Output
    # - Cell_param describing the cell

    # Goes back to the begining of the file if we're not there
    seekstart( file_io )

    # Parse the first line of the file with " " deliminator
    keywords = utils.getLineElements( file_io )

    # Initialize lengths vector
    lengths=zeros(Real,3)

    # Initialize angles vector
    angles=zeros(Real,3)

    # Loop over dimension
    for i=1:3
        # Parse strings into float
        lengths[i] = parse( Float64, keywords[i+2] ) # For lengths for element 3-5
        angles[i]  = parse( Float64, keywords[i+5] ) # For angles for element 6-8
    end

    # Converts lengths and angles into Cell_param and returns it
    return cell_mod.Cell_param( lengths, angles )
end
# Read cell info from .res file using string for the path of the file
function extractCellInfo( file_path::T1 ) where { T1 <: AbstractString }
    # Argument
    # - file_path: path to the .res file
    # Output
    # - cell : contains Cell_param with cell information

    # Opens file
    handle_in = open( file_path )

    # Get cell information
    cell = extractCellInfo( handle_in )

    # Close file
    close( handle_in )

    # Returns Cell_param
    return cell
end
#-------------------------------------------------------------------------------

# Get number of atoms in the structure of the file
#-------------------------------------------------------------------------------*
# Get number of atoms in the file, using IO handler for file
function getNbAtoms( handle_in::T1 ) where { T1 <: IO }
    # Argument
    # - handle_in: IO handler
    # Output
    # - number of atoms in the structure of the file

    # Loop back to the begining to ensure we are there
    seekstart(handle_in)

    # Initialize number of lines
    nb_lines = 0

    # Loop as long as the file can be read
    while ! eof( handle_in )
        # skip line
        readline( handle_in )
        # Increments counter
        nb_lines += 1
    end

    # Number of atoms is number of lines in the file minus 4
    return nb_lines-4
end
# Get number of atoms in the file, using file path for the file
function getNbAtoms( file_path::T1 ) where { T1 <: AbstractString }
    # Argument
    # - file_path: path to the .res file
    # Output
    # Number of atoms in the file
    # OR false if something went wrong with the file

    # Check that the file exists
    if ! isfile( file_path )
        # If the file does not exists, send message and returns false
        print("No .res file at ",file_path,"\n")
        return false
    end

    # Initialize number of lines
    nb_lines=0

    # Opens file
    handle_in = open( file_path )

    # Loop as long as the file can be read
    while ! eof( handle_in )
        # Read line in emptyness
        readline( handle_in )

        # Increments number of lines
        nb_lines += 1
    end

    # Check that the file is not empty or truncated
    if nb_lines - 4 <= 0
        # If it is, sends a message and returns false
        print("File ",file_path," is empty or not valid.\n")
        return false
    end

    # Closing input file
    close(handle_in)

    # Returns number of lines minus 4 for the number of atoms
    return nb_lines-4
end
#-------------------------------------------------------------------------------*


# Extracting atomic informations from file
#-------------------------------------------------------------------------------*
# Read atoms info using IO handler for the file
function extractAtomsInfo( handle_in::T1 ) where { T1 <: IO }
    # Argument
    # - handle_in: IO handler of the target file
    # Output:
    # atoms: AtomList containing atoms information

    # Get number of atoms of the structure
    nb_atoms = getNbAtoms( handle_in )

    # Put read cursor back to the begining of the file
    seekstart( handle_in )

    # Skips 3 first lines
    utils.skipLines( handle_in, 3 )

    # Initialize vectors for positions, names and index
    positions = zeros( Real, nb_atoms, 3 )           # positions
    names = Vector{AbstractString}(undef, nb_atoms)  # names
    index = zeros( Int, nb_atoms )                   # index

    # Loop over atoms
    for atom=1:nb_atoms
        # Parse line with " " deliminator
        keywords = utils.getLineElements( handle_in )

        # The first element is the name of the atom
        names[atom] = keywords[1]

        # Loop over dimensions
        for i=1:3
            # Get elements 2-4 as positions, casts strings into float for positions
            positions[atom,i] = parse(Float64, keywords[ 2+i ] )
        end

        # The loop number gives the index of the atom
        index[atom] = atom
    end

    # Extract the cell information, and converts it into matrix form
    cell = cell_mod.params2Matrix( extractCellInfo( handle_in ) )

    # Transforms the reduced positions into cartesian positions
    positions = cell_mod.getTransformedPosition( positions, cell  )

    # Create AtomList
    atoms = AtomList( names, index, positions )

    # Wraps the atoms if need be
    atoms = cell_mod.wrap( atoms, cell )

    # Returns the AtomList with all informations
    return atoms
end
# Read atoms info using file path for the file
function extractAtomsInfo( file_path::T1 ) where { T1 <: AbstractString }
    # Argument
    # - file_path: path to the target .res file
    # Output
    # - atoms: AtomList with the atomic information

    # Check if file exists
    if ! isfile( file_in )
        # If file does not exists, sends message and return false
        print("No .res file at: ",file_path," !\n")
        return false
    end

    # Opens file
    handle_in = open( file_in)

    # Get atomic information
    atoms = extractAtomsInfo( handle_in )

    # Close file
    close( handle_in )

    # Returns AtomList with atomic information
    return atoms
end
#-------------------------------------------------------------------------------

# Reads *.res file
#-------------------------------------------------------------------------------
function readRes( file_path::T1 ) where { T1 <: AbstractString }
    # Argument:
    # - file_path: path to the input file
    # Output
    # - atoms: AtomList with atomic information
    # - cell: Cell_param with information on the cell

    # Check that the file exists
    if ! isfile( file_path )
        # If it does not exists, sends message and returns false
        print("No .res file at ",file_path," !\n")
        return false
    end

    # Opens file
    handle_in = open( file_path )

    # Get cell information
    cell = extractCellInfo( handle_in )

    # Get Atomic Information
    atoms = extractAtomsInfo( handle_in )

    # Returns AtomList and Cell_parma containing Atomic and cell information on the structure
    return atoms, cell
end
#-------------------------------------------------------------------------------


end
