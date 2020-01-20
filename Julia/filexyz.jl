module filexyz

# filexyz
# Contains functions to deal with xyz file:
# - Functions to get the number of steps in a file
# - Functions to read the file into an AtomList structure (see atom_mod.jl)
# - Functions to write an .xyz file into a given filepath from an AtomList structure

# TODO
# - Support for non AtomList Object
# - Support for AtomMolList
# - Add More write functions? 

export getNbSteps
export readFastFile, readStep, readXYZ
export writeXYZ

using atom_mod
using cell_mod
using cube_mod


# Counts the nb of steps contained in file
#------------------------------------------------------------------------------
function getNbStep( file_path::T1 ) where { T1 <: AbstractString }

  if ! isfile( file_path )
    print("File TRAJEC.xyz does not exists at ",file_path,"\n")
    return false
  end

  nb_lines=0
  nb_atoms=0
  file_in = open( file_path )
  while ! eof(file_in)
    if nb_step == 0
      nb_atoms = parse(Float64,split( readline( file_in ) )[1] )
    else
      temp = readline( file_in )
    end
    nb_lines += 1
  end
  close(file_in)

  if nb_lines % (nb_atoms+2) != 0
    print("File TRAJEC.xyz is corrupted!\n")
    return false
  end

  return Int(nb_lines/(nb_atoms+2))
end
#------------------------------------------------------------------------------

# Reading file by not managing datas
# - requires ability to load both all lines of the file
#   and the array of position of all the trajectory simultaneously
#   fastest implementation, but also more memory hungry
function readFastFile( file::T1 ) where { T1 <: AbstractString }

  # Checking file
  if ! isfile(file)
    return sim=Vector{ atom_mod.AtomList }( undef, 0 ), false
  end

  #--------------
  # Reading file
  #----------------------
  file2=open(file);
  lines=readlines(file2);
  close(file2);
  #------------------------

  #------------------------
  # Basic data about files
  #-----------------------------------------
  if size(lines)[1] == 0
    print("File is empty\n")
    return sim=Vector{ atom_mod.AtomList }( undef, 0 )
  end
  nb_atoms=parse(Int64,split(lines[1])[1])
  if size(lines)[1]/(nb_atoms+2) - trunc(size(lines)[1]/(nb_atoms+2)) > 0.00000000001
    print("ERROR: Problem inside the file! (number of atoms anounced and given do not match).")
    return sim=Vector{ atom_mod.AtomList }( undef, 0 )
  end
  nb_steps=Int(size(lines)[1]/(nb_atoms+2))
  #------------------------------------------

  # Reading file into array
  #------------------------------------------------
  sim=Vector{ atom_mod.AtomList }( undef, nb_steps )
  for step=1:nb_steps
      atom_list = atom_mod.AtomList( nb_atoms )
      for atom=1:nb_atoms
          line_nb=Int((step-1)*(nb_atoms+2)+atom+2)
          line_content=split( lines[line_nb] )
          atom_list.names[atom] = line_content[1]
          atom_list.index[atom] = atom
          for pos=1:3
              atom_list.positions[ atom, pos ] = parse( Float64, line_content[ pos+1 ] )
          end
      end
      sim[step]=atom_list
  end
  #------------------------------------------------

  return sim, true
end

# Read a single step
function readStep( file_hand::T1 ) where { T1 <: IO }
  # Get number of atoms
  line_temp=readline(file_hand)
  nb_atoms = parse(Int,split(line_temp)[1])
  # Atoms
  atoms = AtomList(nb_atoms)
  # Reading comment line
  readline(file_hand)
  # Loop over atoms
  for i=1:nb_atoms
    line=split(readline(file_hand))
    atoms.names[i]=line[1]
    atoms.index[i]=i
    for j=1:3
      atoms.positions[i,j] = parse(Float64,line[j+1])
    end
  end
  return atoms
end

# Advanced the cursor one step forward
function skipStep( file_hand::T1 ) where { T1 <: IO }
  # Get number of atoms
  nb_atoms = parse(Int,split(readline(file_hand))[1])
  # Reading comment line
  readline(file_hand)
  # Loop over atoms
  for i=1:nb_atoms
    line=split(readline(file_hand))
  end
  return
end

# Read a whole file with memory management
function readXYZ( file::T1 ) where { T1 <: AbstractString }
  nb_steps=getNbSteps(file)
  file_hand=open(file)
  atoms_sim=Vector{AtomList}(undef,nb_steps)
  for i=1:nb_steps
    atoms_sim[1] = readStep(file_hand)
  end
  close(file_hand)
  return atoms_sim
end

# Read a file with a given stride
function readXYZ( file::T1, stride::T2 ) where { T1 <: AbstractString, T2 <: Int }
  nb_steps=getNbSteps(file)
  file_hand=open(file)
  atoms_sim=Vector{AtomList}(Int(trunc(nb_steps/stride)))
  j=1
  for i=1:nb_steps
    if i % stride == 0
      atoms_sim[j]=readStep(file_hand)
      j += 1
    else
      skipStep(file_hand)
    end
  end
  return atoms_sim
end

# read a whole file with a stride and starting at a specified steps
function readXYZ( file::T1, stride::T2 , start_step::T3 ) where { T1 <: AbstractString, T2 <: Int , T3 <: Int }
  nb_steps=getNbSteps(file)
  file_hand=open(file)
  atoms_sim=Vector{AtomList}(Int(trunc((nb_steps-start_step)/stride)))
  for i=1:start_step
    readEmpty(file_hand)
  end
  j=1
  for i=start_step:nb_steps-stride
    if i % stride == 0
      atoms_sim[j]=readStep(file_hand)
      j += 1
    else
      readEmpty(file_hand)
    end
  end
  return atoms_sim
end
#--------------------------------------------------------------------------------

# Writing File To Disk
#==============================================================================#
function writeXYZ( file_handle::T1, atoms::T2 ) where { T1 <: IOStream, T2 <: atom_mod.AtomList }
  nb_atoms=size(atoms.names)[1]
  Base.write(file_handle,string(nb_atoms,"\n"))
  Base.write(file_handle,string("STEP: X\n"))
  for i=1:nb_atoms
    Base.write( file_handle, string( atoms.names[i] , " "  ) )
    for j=1:3
      Base.write( file_handle, string( atoms.positions[i,j] , " " ) )
    end
    Base.write( file_handle, "\n" )
  end
  return true
end
function writeXYZ( file_handle::T1, traj::Vector{T2} ) where { T1 <: IOStream, T2 <: atom_mod.AtomList }
  nb_step=size(traj)[1]
  for step=1:nb_step
    temp=writeXYZ( file_handle, traj[step] )
  end
  return true
end
function writeXYZ( file_path::T1, atoms::T2 ) where { T1 <: AbstractString, T2 <: atom_mod.AtomList }
  file_out = open( file_path, "w" )
  nb_atoms=size(atoms.names)[1]
  Base.write( out, string( nb_atoms, "\n"))
  Base.write( out, string( "STEP: X\n" ) )
  for i=1:nb_atoms
    write( file_out, string( atoms.names[i], " " ) )
    for j=1:3
      write( file_out, string( atoms.positions[i,j] , " ") )
    end
    write( file_out, "\n" )
  end
  close( file_out )
  return true
end
function writeXYZ( file::T1, traj::Vector{T2} ) where { T1 <: AbstractString, T2 <: atom_mod.AtomList }
  nb_step=size(traj)[1]
  file_out = open(file,"w")
  for step=1:nb_step
    writeXYZ( file_out, traj[step] )
  end
  close(file_out)
  return true
end
#==============================================================================#

end
