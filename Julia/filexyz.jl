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
function getNbStepAtoms( file_path::T1 ) where { T1 <: AbstractString }

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

  return Int(nb_lines/(nb_atoms+2)), nb_atoms
end
#------------------------------------------------------------------------------


# Reading XYZ files
#==============================================================================#
function readFileAtomList( file_path::T1 ) where { T1 <: AbstractString }

  #-----------------------------------------------
  nb_step = getNbSteps( file_path )
  if nb_step == false
    return false
  end
  #-----------------------------------------------

  # Init output structure
  #------------------------------------------------
  traj=Vector{ atom_mod.AtomList }( undef, nb_steps )
  #-------------------------------------------------

  # Reading
  #------------------------------------------------------
  file_in = open( file_path )
  for step=1:nb_steps
      nb_atoms=parse( Int64, split( readline( file_in ) )[1] )
      temp=readline( file_in )
      traj[step] = atom_mod.AtomList( nb_atoms )
      for atom=1:nb_atoms
          keywords=split( readline(file_in) )
          traj[step].names[atom] = line_keywords[1]
          traj[step].index[atom] = atom
          for i=1:3
              traj[step].positions[ atom, pos ] = parse( Float64, keywords[ i+1 ] )
          end
      end
  end
  #------------------------------------------------

  return traj
end
function readFileAtomList( file_path::T1, stride_::T2 ) where { T1 <: AbstractString, T2 <: Int }

  #-----------------------------------------------
  nb_step = getNbSteps( file_path )
  if nb_step == false
    return false
  end
  #-----------------------------------------------

  # Init output structure
  #------------------------------------------------
  traj=Vector{ atom_mod.AtomList }( undef, nb_steps )
  #-------------------------------------------------

  # Reading
  #------------------------------------------------------
  file_in = open( file_path )
  for step=1:nb_steps
      nb_atoms=parse( Int64, split( readline( file_in ) )[1] )
      temp=readline( file_in )
      traj[step] = atom_mod.AtomList( nb_atoms )
      for atom=1:nb_atoms
          keywords=split( readline(file_in) )
          traj[step].names[atom] = line_keywords[1]
          traj[step].index[atom] = atom
          for i=1:3
              traj[step].positions[ atom, pos ] = parse( Float64, keywords[ i+1 ] )
          end
      end
  end
  #------------------------------------------------

  return traj
end
#==============================================================================#

# Writing XYZ File To Disk (currently solely from AtomList structures )
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
