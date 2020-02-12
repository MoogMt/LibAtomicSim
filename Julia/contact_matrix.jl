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

# Getting bonds
#-------------------------------------------------------------------------------
function printBonds( atoms::T1, cell::T2, index::T3 , cut_off::T4 ) where { T1 <: atom_mod.AtomList, T2 <: cell_mod.Cell_param , T3 <: Int , T4 <: Real }
  nb_atoms=size(atoms.names)[1]
  for i=1:nb_atoms
    dist=cell_mod.distance(atoms,cell,index,i)
    if dist < cut_off
      print("Bond with ", i, "\n");
    end
  end
  return
end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
function readMatrix( input::T1 , nb_atoms::T2) where { T1 <: IO, T2 <: Int }
  matrix=zeros(nb_atoms,nb_atoms)
  for i=1:nb_atoms
    line=split(readline(input))
    for j=1:nb_atoms
      matrix[i,j] = parse( Float64, line[j] )
    end
  end
  if split(readline(input))[1] == "END"
    return matrix
  else
    print("Problem while reading file")
    return
  end
end
function readMatrix( file::T1 ) where { T1 <: AbstractString }
  file=open(file)
  nb_steps=Int(split(readline(file))[1])
  close(file)
  return nb_steps
end
function readMatrix( file::T1, step::T2 ) where { T1<: AbstractString, T2<: Int }
  file=open(file)
  line=split(readline(file))

  nb_steps=parse(Int,line[1]);
  nb_atoms=parse(Int,line[2]);
  if nb_steps < step
    return false
  end
  for i=1:step
    for j=1:nb_atoms
      line=readline(file)
    end
  end
  matrix=zeros(nb_atoms,nb_atoms)
  for i=1:nb_atoms
    line=split(readline(file))
    for j=1:nb_atoms
      matrix[i,j]=parse(Float64,line[j])
    end
  end
  close(file)
  return matrix
end
#-------------------------------------------------------------------------------

#--------------------------
function writeMatrix( file::T1, matrix::Array{T2} ) where { T1 <: IO , T2 <: Real }
    nb_atoms=size(matrix)[1]
    for i=1:nb_atoms
        for j=1:nb_atoms
            write(file, string(matrix[i,j]," ") )
        end
        write(file,"\n")
    end
    return
end
#--------------------------

end
