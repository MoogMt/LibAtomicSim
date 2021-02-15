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
    handle_in = open( file_path )

    # Loop over lines as long as there are lines to read
    while ! eof(handle_in)
        line=readline( handle_in )
        # The first line of the file gives the number of atoms
        if nb_lines == 0
            # Parse the string with number of atoms from string to float
            nb_atoms = parse(Float64, split( line )[1] )
        end
        # Increments counter of lines
        nb_lines += 1
    end

    # Close file
    close(handle_in)

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
    handle_in = open( file_path )

    # Loop over line as long as possible
    while ! eof(handle_in)
        # Read line
        line = readline( handle_in )

        # The first line contains the number of atoms
        if nb_step == 0
            # Parsing the string into a float for number of atoms
            nb_atoms = parse(Float64, split( line )[1] )
        end

        # Increments line counter
        nb_lines += 1
    end

    # Closes file
    close(handle_in)

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
    handle_in = open( file_path )

    # Get the number of atoms in the structure as the first element of the first line
    # + cast String into Int
    nb_atoms = parse(Int64, split( readline( handle_in ) )[1] )

    # Skipping comment line
    temp=readline( handle_in )

    # Initialize AtomList for output
    atoms = atom_mod.AtomList( nb_atoms )

    # Loop over atoms
    for atom=1:nb_atoms
        # Reads line and split with " " deliminator
        keywords=split( readline(handle_in) )

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
        return false
    end

    # If number of step is 1, use appropriate function
    if nb_step == 1
        return readStructureAtomList( file_path )
    end

    # Initialize vector of AtomList for trajectory
    traj = Vector{ atom_mod.AtomList }(undef, nb_step )

    # Opens input file
    handle_in = open( file_path )

    # Loop over steps
    for step=1:nb_step
        # Get number of atoms as first element of first line
        nb_atoms = parse( Int64, split( readline( hanle_in ) )[1] )

        # Skip line command
        temp = readline( handle_in )

        # Initialize
        traj[step] = atom_mod.AtomList( nb_atoms )

        # Loop over atoms
        for atom=1:nb_atoms
            # Read line and parse with " "
            keywords = split( readline( handle_in ) )

            # First element of the line is the Atom
            traj[step].names[atom] = keywords[1]

            # Index is the number of atom
            traj[step].index[atom] = atom

            # Loop over dimensions
            for i=1:3
                # Parse the position as elements 2-4
                traj[step].positions[ atom, i ] = parse( Float64, keywords[ i+1 ] )
            end
        end
    end

    # Closing line
    close( handle_in )

    # Return vector of AtomList
    return traj
end
# Reads *.xyz file into vector of AtomList for trajectory using a stride
function readFileAtomList( file_path::T1, stride_::T2 ) where { T1 <: AbstractString, T2 <: Int }
    # Argument
    # - file_path: path to the *.xyz file
    # - stride_ : (int) stride with which to read the trajectory
    # Output
    # - traj: vectors of AtomList for the traj
    # OR false

    # Get the initial number of steps in the file
    nb_step_origin = getNbSteps( file_path )
    # Check that the file exists and that the file is not empty
    if nb_step_origin == false || nb_step_origin == 0
        # Sends a message and return false if problem
        print("File is empty or does not exists.\n")
        return false
    end

    # Compute number of steps post striding
    nb_step = utils.nbStepStriding( nb_step_origin, stride_ )

    # Initialize vector of AtomList for trajectory
    traj = Vector{ atom_mod.AtomList }( undef, nb_step )

    # Opens file
    handle_in = open( file_path )

    # Initialize step counter
    count_step=1

    # Loop over all steps
    for step=1:nb_step_origin
        # Get the number of atoms as the first element of the frame
        nb_atoms=parse( Int64, split( readline( handle_in ) )[1] )

        # Skip comment line
        temp = readline( handle_in )

        # Check whether or not stride
        if (step-1) % stride_ == 0
            # Initialize AtomList current step
            traj[count_step] = atom_mod.AtomList( nb_atoms )

            # Loop over atoms
            for atom=1:nb_atoms
                # Read and parse line with " " deliminator
                keywords = split( readline( handle_in ) )

                # Get the name of the atom as first element of the line
                traj[count_step].names[atom] = keywords[1]

                # Get the index number as the atom number
                traj[count_step].index[atom] = atom

                # Loop over dimensions
                for i=1:3
                    # Parse elements 2-4 as atomic positions
                    traj[count_step].positions[ atom, i ] = parse( Float64, keywords[ i+1 ] )
                end
            end

            # Add steps
            count_step += 1
        else
            # If we are in strided frame we skip all atomic positions
            for i=1:nb_atoms
                # Skipping atom lines
                temp = readline( handle_in )
            end
        end
    end

    # Close input file
    close( handle_in )

    # Returns vector of AtomList
    return traj
end
# Reads *.xyz file into vector of AtomList for trajectory
# -> using a stride
# -> an initial number of step ignored
function readFileAtomList( file_path::T1, stride_::T2, nb_ignored::T3 ) where { T1 <: AbstractString, T2 <: Int, T3 <: Int }
    # Argument
    # - file_path: path to the *.xyz file
    # - stride_ : (int) stride with which to read the trajectory
    # - nb_ignored: (int) number of step to ignore at the begining
    # Output
    # - traj: vectors of AtomList for the traj
    # OR false if there is problem with the file

    # Get the initial number of steps in the file
    nb_step_origin = getNbSteps( file_path )
    # Check that the file exists and that the file is not empty
    if nb_step_origin == false || nb_step_origin == 0
        # Sends a message and return false if problem
        print("File is empty or does not exists.\n")
        return false
    end

    # Get number of steps with striding
    nb_step = utils.nbStepStriding( nb_step_origin - nb_ignored, stride_ )

    # Vector of AtomList, for the trajectory
    traj = Vector{ atom_mod.AtomList }( undef, nb_step )

    # Opens file
    handle_in = open( file_path )

    # Initialize step counter
    count_step=1

    # Loop over steps
    for step=1:( nb_step_origin - nb_ignored )
        # We care about the number of atoms to skip...
        nb_atoms = parse( Int64, split( readline( handle_in ) )[1] )

        # Skipping comment line
        temp = readline( handle_in )

        # Check if step is strided or not
        if ( step - 1 ) % stride_ == 0
            # Initialize AtomList
            traj[count_step] = atom_mod.AtomList( nb_atoms )

            # Loop over atoms
            for atom=1:nb_atoms
                # Reading line and parsing with " "
                keywords=split( readline( handle_in ) )

                # First element of the line is atom name
                traj[count_step].names[atom] = keywords[1]

                # Putting index atom
                traj[count_step].index[atom] = atom

                # Loop over dimensions
                for i=1:3
                    # Elements 2-4 are atomic positions, parsing string into floats
                    traj[count_step].positions[ atom, i ] = parse( Float64, keywords[ i+1 ] )
                end
            end
            # Increment step counter
            count_step += 1
        else
            # Strided step, skipping
            for atom=1:nb_atoms
                # Skip all atomic lines
                temp = readline( handle_in )
            end
        end
    end

    # Close file
    close( handle_in )

    # Returns a vector of AtomList
    return traj
end
# Reads *.xyz file into vector of AtomList for trajectory
# -> using a stride
# -> an initial number of step ignored
# -> a number of maximum step to count
function readFileAtomList( file_path::T1, stride_::T2, nb_ignored::T3, nb_max::T4 ) where { T1 <: AbstractString, T2 <: Int, T3 <: Int, T4 <: Int }
    # Argument
    # - file_path: path to the *.xyz file
    # - stride_ : (int) stride with which to read the trajectory
    # - nb_ignored: (int) number of step to ignore at the begining
    # - nb_max: (int) max number of steps
    # Output
    # - traj: vectors of AtomList for the traj
    # OR false if there is problem with the file

    # Get the initial number of steps in the file
    nb_step_origin = getNbSteps( file_path )
    # Check that the file exists and that the file is not empty
    if nb_step_origin == false || nb_step_origin == 0
        # Sends a message and return false if problem
        print("File is empty or does not exists.\n")
        return false
    end

    # Compute number of steps post striding
    nb_step = utils.nbStepStriding( nb_step_origin - nb_ignored, stride_ )

    # If nb_max is negative, returns false
    if nb_max <= 0
        print("nb_max must be positive!\n")
        return false
    end

    # Initialize vector of AtomList for the trajectory
    traj=Vector{ atom_mod.AtomList }( undef, nb_max )

    # Opens file
    handle_in = open( file_path )

    # Initialize step counter
    count_step=1

    # Loop over steps
    for step=1:( nb_step_origin - nb_ignored )
        # We get the number of atoms as first element of the first line of frame
        nb_atoms = parse(Int64, split( readline( handle_in ) )[1] )

        # Skip comment line
        temp = readline( handle_in )

        # Check whether or not to stride step
        if ( step - 1 ) % stride_ == 0
            # Initialize current step AtomList
            traj[count_step] = atom_mod.AtomList( nb_atoms )

            # Loop over atoms
            for atom=1:nb_atoms
                # Reads atom line with " " deliminator
                keywords = split( readline( handle_in ) )

                # The first element of the line is the atom name
                traj[count_step].names[atom] = keywords[1]

                # The index of the atom is the number of the loop
                traj[count_step].index[atom] = atom

                # Loop over dimension
                for i=1:3
                    # The elements 2-4 are casted from String to Float as atomic positions
                    traj[count_step].positions[ atom, i ] = parse( Float64, keywords[ i+1 ] )
                end
            end

            # If the step counter is larger than nb_max, stopping the loop
            if count_step >= nb_max
                break
            end

            # Increments step counter
            count_step += 1
        else
            # Not counting skipped steps
            for atom=1:nb_atoms
                # Skipping atomic line
                temp = readline( handle_in )
            end
        end
    end

    # Close file
    close( handle_in )

    # Returns Vector of AtomList for the trajectory
    return traj
end
#------------------------------------------------------------------------------

# Writing XYZ File To Disk (currently solely from AtomList structures )
#------------------------------------------------------------------------------
# Write single structure *.xyz using AtomList using IO handler for file
function writeXYZ( file_handle::T1, atoms::T2 ) where { T1 <: IOStream, T2 <: atom_mod.AtomList }
    # Argument
    # - file_handle: IO handler for the file
    # - atoms: AtomList with atomic positions
    # Output
    # - Bool: whether the file was written successfully

    # Get number of atoms
    nb_atoms = size( atoms.names )[1]

    # Write atom number line
    Base.write(file_handle,string(nb_atoms,"\n"))

    # Write comment line
    Base.write(file_handle,string("STEP: X\n"))

    # Loop over atoms
    for atom=1:nb_atoms
        # Write atom name
        Base.write( file_handle, string( atoms.names[i] , " "  ) )

        # Loop over dimensions
        for j=1:3
            Base.write( file_handle, string( atoms.positions[i,j] , " " ) )
        end

        # Write end line
        Base.write( file_handle, "\n" )
    end

    # Return True if it was right
    return true
end
# Write trajectory *.xyz using vector of AtomList using IO handler for file
function writeXYZ( file_handle::T1, traj::Vector{T2} ) where { T1 <: IOStream, T2 <: atom_mod.AtomList }
    # Argument
    # - file_handle: IO handler for the file
    # - traj: vector of AtomList with atomic positions
    # Output
    # - Bool: whether the file was written successfully

    # Get number of steps
    nb_step=size(traj)[1]

    # Loop over steps
    for step=1:nb_step
        # Write steps
        temp = writeXYZ( file_handle, traj[step] )

        # If something went wrong, returns false
        if ! temp
            return false
        end
    end

    # Returns True if all went well
    return true
end
# Write single structure *.xyz using vector of AtomList using file_path
function writeXYZ( file_path::T1, atoms::T2 ) where { T1 <: AbstractString, T2 <: atom_mod.AtomList }
    # Argument
    # - file_path: string with the file path of the output file
    # - atoms: AtomList with atomic positions and names
    # Output
    # - Bool: whether or not to write file

    # Opens file
    file_out = open( file_path, "w" )

    # Get number of atoms of the structure
    nb_atoms = size( atoms.names )[1]

    # Write number of atoms
    write( file_out, string( nb_atoms, "\n"))

    # Write comment line
    write( file_out, string( "STEP: X\n" ) )

    # Loop over atoms
    for atom=1:nb_atoms
        # Write atom name at the begining of the line
        write( file_out, string( atoms.names[atom], " " ) )

        # Loop over dimensions
        for i=1:3
            # Write atomic positions
            write( file_out, string( atoms.positions[atom,i] , " ") )
        end

        # Write end of line
        write( file_out, "\n" )
    end

    # Closes file
    close( file_out )

    # Returns true if all went well
    return true
end
# Write trajectory *.xyz using vector of AtomList using file_path
function writeXYZ( file_path::T1, traj::Vector{T2} ) where { T1 <: AbstractString, T2 <: atom_mod.AtomList }
    # Argument
    # - file_path: path of the file to write data in
    # - traj: vector of AtomList
    # Output:
    # - Bool: whether or not writting was succesful

    # Getting number of step
    nb_step = size( traj )[1]

    # Opening output file
    handle_out = open(file,"w")

    # Loop over steps
    for step = 1:nb_step
        # Writing AtomList to file
        writeXYZ( handle_out, traj[step] )
    end

    # Close file
    close( handle_out )

    # Returns true if something went ok
    return true
end
#------------------------------------------------------------------------------

end
