module contact_matrix

export buildMatrix, readMatrix
export getBonded, computeMatrix, writeMatrix

using atom_mod
using cell_mod
using graph

# Building Matrix
#-------------------------------------------------------------------------------
function buildMatrix( atoms::T1, cell::T2 ) where { T1 <: atom_mod.AtomList , T2 <: cell_mod.Cell_param }
    nb_atoms=size(atoms.names)[1]
    matrix=zeros(nb_atoms,nb_atoms)
    for i=1:nb_atoms
        for j=1:nb_atoms
            dist=cell_mod.distance(atoms,cell,i,j)
            matrix[i,j]=dist
            matrix[j,i]=dist
        end
    end
    return matrix
end
function buildMatrix( atoms::T1 , cell::T2, cut_off::T3 ) where { T1 <: atom_mod.AtomList , T2 <: cell_mod.Cell_param, T3 <: Real }
    nb_atoms=size(atoms.names)[1]
    matrix=zeros(Int,nb_atoms,nb_atoms)
    for i=1:nb_atoms
        for j=i+1:nb_atoms
            if cell_mod.distance(atoms,cell,i,j) <= cut_off
                matrix[i,j]=1
                matrix[j,i]=1
            end
        end
    end
    return matrix
end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
function readStepMatrix( handle_in::T1 , nb_atoms::T2) where { T1 <: IO, T2 <: Int }
    matrix=zeros(nb_atoms,nb_atoms)
    for i=1:nb_atoms
        keywords = split( readline( handle_in ) )
        for j=1:nb_atoms
            matrix[i,j] = parse( Float64, keywords[j] )
        end
    end
    return matrix
end
function readMatrix( file::T1, target_step::T2 ) where { T1 <: AbstractString, T2 <: Int }
    handle_in = open( file )
    # Getting info
    keywords = split( readline( file ) )
    nb_step  = parse(Int, keywords[1] )
    nb_atoms = parse(Int, keywords[2] )
    if step > nb_step
        print("The trajectory is only ", nb_step, " long, you're asking for step ", target_step,".\n")
        print("Stopping now!\n")
        return false
    end
    # Reading
    matrix=zeros(nb_atoms,nb_atoms)
    for step=1:nb_step
        if step == target_step
            matrix[ :, : ]  = readStepMatrix( handle_in, nb_atoms )
        else
            readStepMatrix( handle_in, nb_atoms )
        end
    end
    close(file)
    return matrix
end
function readMatrix( file::T1 ) where { T1<: AbstractString }
    handle_in = open( file )
    # Getting info
    keywords = split( readline( file ) )
    nb_step  = parse(Int, keywords[1] )
    nb_atoms = parse(Int, keywords[2] )
    # Reading
    matrix=zeros(nb_step,nb_atoms,nb_atoms)
    for step=1:nb_step
        matrix[ step, :, : ]  = readStepMatrix( handle_in, nb_atoms )
    end
    close(file)
    return matrix
end
#-------------------------------------------------------------------------------

# Writting
#-------------------------------------------------------------------------------
function writeStepMatrix( handle_out::T1, matrix::Array{T2,2} ) where { T1 <: IO , T2 <: Real }
    nb_atoms=size(matrix)[1]
    for atom1=1:nb_atoms
        for atom2=1:nb_atoms
            write( handle_out, string( matrix[ atom1, atom2 ], " " ) )
        end
        write( handle_out , "\n")
    end
    return true
end
function writeMatrix( file_out::T1, matrix::Array{T2,3} ) where { T1 <: AbstractString, T2 <: Real }
    handle_out = open( file_out, "w" )
    nb_step = size( matrix )[1]
    nb_atoms = size( matrix )[2]
    Base.write( handle_out, string( nb_step, " ", nb_atoms, "\n" ) )
    for step=1:nb_atoms
        writeStepMatrix( handle_out, matrix[step,:,:] )
    end
    close( handle_out )
    return true
end
#-------------------------------------------------------------------------------


end
