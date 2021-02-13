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
    handle_in = open( file_path )
    cell = extractCellInfo( handle_in )
    close( handle_in )
    return cell
end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------*
function getNbAtoms( handle_in::T1 ) where { T1 <: IO }
    seekstart(handle_in)
    nb_lines=0
    while ! eof( handle_in )
        test = readline( handle_in )
        nb_lines += 1
    end
    return nb_lines-4
end
function getNbAtoms( file_path::T1 ) where { T1 <: AbstractString }

    #-------------------------------------------
    if ! isfile( file_path )
        print("No .res file at ",file_path,"\n")
        return false
    end
    #-------------------------------------------

    #-------------------------------------------
    nb_lines=0
    handle_in = open( file_path )
    while ! eof( handle_in )
        test = readline( handle_in )
        nb_lines += 1
    end
    close(handle_in)
    #-------------------------------------------

    return nb_lines-4
end
function extractAtomsInfo( handle_in::T1 ) where { T1 <: IO }
    nb_atoms = getNbAtoms( handle_in )
    seekstart( handle_in )
    utils.skipLines( handle_in, 3 )
    positions=zeros( Real, nb_atoms, 3 )
    names = Vector{AbstractString}(undef, nb_atoms)
    index = zeros( Int, nb_atoms )
    for atom=1:nb_atoms
        keywords = utils.getLineElements( handle_in )
        names[atom] = keywords[1]
        for i=1:3
            positions[atom,i] = parse(Float64, keywords[ 2+i ] )
        end
        index[atom] = atom
    end
    cell = cell_mod.params2Matrix( extractCellInfo( handle_in ) )
    positions=cell_mod.getTransformedPosition( positions, cell.matrix  )
    atoms=AtomList( names, index, positions )
    atoms=cell_mod.wrap( atoms, cell )
    return atoms
end
function extractAtomsInfo( file_path::T1 ) where { T1 <: AbstractString }

    #------------------------------------------------
    if ! isfile( file_in )
        print("No .res file at: ",file_path," !\n")
        return false
    end
    #-------------------------------------------------

    #------------------------------------
    handle_in = open( file_in)
    atoms = extractAtomsInfo( handle_in )
    close( handle_in )
    #-------------------------------------

    return atoms
end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
function readRes( file_path::T1 ) where { T1 <: AbstractString }

    if ! isfile( file_path )
        print("No .res file at ",file_path," !\n")
        return false
    end

    handle_in = open( file_path )
    cell = extractCellInfo( handle_in )
    atoms = extractAtomsInfo( handle_in )

    return atoms, cell
end
#-------------------------------------------------------------------------------


end
