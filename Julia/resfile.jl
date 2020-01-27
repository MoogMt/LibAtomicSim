module resfile

using utils
using conversion
using atom_mod
using cell_mod
using filexyz
using pdb

using LinearAlgebra

#-------------------------------------------------------------------------------
function extractCellInfo( file_io::T1 ) where { T1 <: IO}
    # incase not at start
    seekstart( file_io )
    keywords = utils.getLineElements( file_io )
    lengths=zeros(Real,3)
    for i=1:3
        lengths[i] = parse( Float64, keywords[i+2] )
    end
    angles=zeros(Real,3)
    for i=1:3
        angles[i] = parse( Float64, keywords[i+5] )
    end
    return cell_mod.Cell_param( lengths, angles )
end
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
    positions=cell_mod.reduced2Cartesian( positions, cell_mod.params2Matrix( extractCellInfo( handle_in ) ) )
    return AtomList( names, index, positions )
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
