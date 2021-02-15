module filexyz

# Loading necessary modules from LibAtomicSim
using atom_mod
using cell_mod
using cube_mod
using utils

# Description
# Contains functions to deal with xyz file:
# - Functions to get the number of steps in a file
# - Functions to read the file into an AtomList structure (see atom_mod.jl)
# - Functions to write an .xyz file into a given filepath from an AtomList structure

# Exporting functions
export getNbSteps
export readFastFile, readStep, readXYZ
export writeXYZ

# TODO
# - Support for non AtomList Object
# - Support for AtomMolList
# - Add More write functions?
# - Check that all functions are exported

# Getting number of atoms and steps of a *.xyz file
#------------------------------------------------------------------------------
# Get number of steps of a file
function getNbSteps( file_path::T1 ) where { T1 <: AbstractString }
  # Argument:
  # file_path: path to the *.xyz file
  # Output:
  # Int, number of steps in the trajectory

  #
  # Check that the file exists
  if ! isfile( file_path )
    # If file is empty, send a message and returns false
    print("File TRAJEC.xyz does not exists at ",file_path,"\n")
    return false
  end

  # Initialize number of atoms and lines
  nb_lines=0
  nb_atoms=0

  # Open input file
  file_in = open( file_path )

  # Loop over lines as long as there are lines to read
  while ! eof(file_in)
    line=readline( file_in )
    # The first line of the file gives the number of atoms
    if nb_lines == 0
      # Parse the string with number of atoms from string to float
      nb_atoms = parse(Float64, split( line )[1] )
    end
    # Increments counter of lines
    nb_lines += 1
  end

  # Close file
  close(file_in)

  # Check that the file format is ok
  # - There should be nb_atoms+2 lines per step
  if nb_lines % ( nb_atoms + 2 ) != 0
    # If the file format is wrong, sends a message and return false
    print("File TRAJEC.xyz at ",file_path," is corrupted!\n")
    return false
  end

  # Returns number of atoms
  return Int( nb_lines/(nb_atoms + 2 ) )
end
# Gets the number of steps and number of atoms of a file
# -> Assumes that it is a traj and the number of atoms does not change
function getNbStepAtoms( file_path::T1 ) where { T1 <: AbstractString }
  # Argument
  # - file_path: path to the input *.xyz file
  # Output
  # - Int, number of steps
  # - nb_atoms: nb of atoms, Int

  # Check that the input file exists
  if ! isfile( file_path )
    # If not, sends a message and returns false
    print("File TRAJEC.xyz does not exists at ",file_path,"\n")
    return false
  end

  # Initialize number of atoms and lines
  nb_lines=0
  nb_atoms=0

  # Opens the file
  file_in = open( file_path )

  # Loop over line as long as possible
  while ! eof(file_in)
    # Read line
    line = readline( file_in )

    # The first line contains the number of atoms
    if nb_step == 0
      # Parsing the string into a float for number of atoms
      nb_atoms = parse(Float64, split( line )[1] )
    end

    # Increments line counter
    nb_lines += 1
  end

  # Closes file
  close(file_in)

  # Check the file format
  if nb_lines % (nb_atoms+2) != 0
    # If the file format is problematic, sends a message and return false
    print("File TRAJEC.xyz at ",file_path," is probably corrupted!\n")
    return false
  end

  # Returns number of steps and number of atoms in the traj
  return Int(nb_lines/(nb_atoms+2)), nb_atoms
end
#------------------------------------------------------------------------------

# Reading XYZ files
#------------------------------------------------------------------------------
# Reads a single structure from *.xyz file into AtomList
function readStructureAtomList( file_path::T1 ) where { T1 <: AbstractString }
  # Argument
  # - file_path: path to the *.xyz file
  # Output
  # - atoms: AtomList containing the structure
  # OR false if something went wrong

  # Check if the file exists
  if ! isfile(file_path)
    # If the file does not exists, sends a message and return false
    print("No file found at: ",file_path,"\n")
    return false
  end

  # Opens file
  file_in = open( file_path )

  # Get the number of atoms in the structure as the first element of the first line
  # + cast String into Int
  nb_atoms = parse(Int64, split( readline( file_in ) )[1] )

  # Skipping comment line
  temp=readline( file_in )

  # Initialize AtomList for output
  atoms = atom_mod.AtomList( nb_atoms )

  # Loop over atoms
  for atom=1:nb_atoms
    # Reads line and split with " " deliminator
    keywords=split( readline(file_in) )

    # First element is the atom name
    atoms.names[atom] = keywords[1]

    # The atom index is the number of atom
    atoms.index[atom] = atom

    # Loop over dimensions
    for i=1:3
      # The element 2-4 are atomic positions, cast them into float and put them in AtomList positions
      atoms.positions[ atom, i ] = parse( Float64, keywords[ i+1 ] )
    end
  end

  # Return AtomList with structure
  return atoms
end
# Reads *.xyz file into vector of AtomList for trajectory
function readFileAtomList( file_path::T1 ) where { T1 <: AbstractString }
  # Argument
  # - file_path: path to the *.xyz file
  # Output
  # - traj: vector of AtomList, contains the atomic trajectory
  # OR false if something went wrong with the file

  # Get the number of steps in the trajectory
  nb_step = getNbSteps( file_path )
  # If the file is empty or does not exists, sends a message and returns false
  if nb_step == false || nb_step == 0
    # I
    return false
  end

  #-----------------------------------------------
  if nb_step == 1
    return readStructureAtomList( file_path )
  end
  #-----------------------------------------------

  # Init output structure
  #------------------------------------------------
  traj=Vector{ atom_mod.AtomList }( undef, nb_step )
  #-------------------------------------------------

  # Reading
  #------------------------------------------------------
  file_in = open( file_path )
  for step=1:nb_step
      nb_atoms=parse( Int64, split( readline( file_in ) )[1] )
      temp=readline( file_in )
      traj[step] = atom_mod.AtomList( nb_atoms )
      for atom=1:nb_atoms
          keywords=split( readline(file_in) )
          traj[step].names[atom] = keywords[1]
          traj[step].index[atom] = atom
          for i=1:3
              traj[step].positions[ atom, i ] = parse( Float64, keywords[ i+1 ] )
          end
      end
  end
  #------------------------------------------------

  return traj
end
function readFileAtomList( file_path::T1, stride_::T2 ) where { T1 <: AbstractString, T2 <: Int }

  #-----------------------------------------------
  nb_step_origin = getNbSteps( file_path )
  if nb_step_origin == false
    return false
  end
  #-----------------------------------------------

  # Init output structure
  #------------------------------------------------
  nb_step = utils.nbStepStriding(nb_step_origin, stride_ )
  traj    = Vector{ atom_mod.AtomList }( undef, nb_step )
  #-------------------------------------------------

  # Reading
  #------------------------------------------------------
  file_in = open( file_path )
  count_step=1
  for step=1:nb_step_origin
      # We care about the number of atoms to skip...
      nb_atoms=parse( Int64, split( readline( file_in ) )[1] )
      # This line is always skipped
      temp=readline( file_in )
      if (step-1) % stride_ == 0
          traj[count_step] = atom_mod.AtomList( nb_atoms )
          for atom=1:nb_atoms
              keywords=split( readline(file_in) )
              traj[count_step].names[atom] = keywords[1]
              traj[count_step].index[atom] = atom
              for i=1:3
                  traj[count_step].positions[ atom, i ] = parse( Float64, keywords[ i+1 ] )
              end
          end
          count_step += 1
      else
          # Not counting those steps
          for i=1:nb_atoms
              temp=readline( file_in )
          end
      end
  end
  #------------------------------------------------

  return traj
end
function readFileAtomList( file_path::T1, stride_::T2, nb_ignored::T3 ) where { T1 <: AbstractString, T2 <: Int, T3 <: Int }

  #-----------------------------------------------
  nb_step_origin = getNbSteps( file_path )
  if nb_step_origin == false
    return false
  end
  #-----------------------------------------------

  # Init output structure
  #------------------------------------------------
  nb_step = utils.nbStepStriding( nb_step_origin-nb_ignored, stride_ )
  traj    = Vector{ atom_mod.AtomList }( undef, nb_step )
  #-------------------------------------------------

  # Reading
  #------------------------------------------------------
  file_in = open( file_path )
  count_step=1
  for step=1:nb_step_origin-nb_ignored
      # We care about the number of atoms to skip...
      nb_atoms=parse( Int64, split( readline( file_in ) )[1] )
      # This line is always skipped
      temp=readline( file_in )
      if (step-1) % stride_ == 0
          traj[count_step] = atom_mod.AtomList( nb_atoms )
          for atom=1:nb_atoms
              keywords=split( readline(file_in) )
              traj[count_step].names[atom] = keywords[1]
              traj[count_step].index[atom] = atom
              for i=1:3
                  traj[count_step].positions[ atom, i ] = parse( Float64, keywords[ i+1 ] )
              end
          end
          count_step += 1
      else
          # Not counting those steps
          for i=1:nb_atoms
              temp=readline( file_in )
          end
      end
  end
  #------------------------------------------------

  return traj
end
function readFileAtomList( file_path::T1, stride_::T2, nb_ignored::T3, nb_max::T4 ) where { T1 <: AbstractString, T2 <: Int, T3 <: Int, T4 <: Int }

  #-----------------------------------------------
  nb_step_origin = getNbSteps( file_path )
  if nb_step_origin == false
    return false
  end
  #-----------------------------------------------

  # Init output structure
  #------------------------------------------------
  nb_step = utils.nbStepStriding( nb_step_origin-nb_ignored, stride_ )
  if nb_max > nb_step
      print("nb_max is too large, maximum value is ",nb_step,"\n")
  end
  if nb_max <= 0
      print("nb_max must be positive!\n")
  end
  traj=Vector{ atom_mod.AtomList }( undef, nb_max )
  #-------------------------------------------------

  # Reading
  #------------------------------------------------------
  file_in = open( file_path )
  count_step=1
  for step=1:nb_step_origin-nb_ignored
      # We care about the number of atoms to skip...
      nb_atoms=parse( Int64, split( readline( file_in ) )[1] )
      # This line is always skipped
      temp=readline( file_in )
      if (step-1) % stride_ == 0
          traj[count_step] = atom_mod.AtomList( nb_atoms )
          for atom=1:nb_atoms
              keywords=split( readline(file_in) )
              traj[count_step].names[atom] = keywords[1]
              traj[count_step].index[atom] = atom
              for i=1:3
                  traj[count_step].positions[ atom, i ] = parse( Float64, keywords[ i+1 ] )
              end
          end
          if count_step >= nb_max
              break
          end
          count_step += 1
      else
          # Not counting those steps
          for i=1:nb_atoms
              temp=readline( file_in )
          end
      end
  end
  #------------------------------------------------

  return traj
end
#------------------------------------------------------------------------------

# Writing XYZ File To Disk (currently solely from AtomList structures )
#------------------------------------------------------------------------------
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
  write( file_out, string( nb_atoms, "\n"))
  write( file_out, string( "STEP: X\n" ) )
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
#------------------------------------------------------------------------------

end
