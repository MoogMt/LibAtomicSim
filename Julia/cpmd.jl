module cpmd

# Load necessary modules from LibAtomicSim
using conversion
using utils
using filexyz
using press_stress
using atom_mod

# Description
# Set of fucntions used to deal with CPMD specific input and output files

# Exporting function
export readInputTimestep, readIntputStrideStress, readIntputStrideTraj
export readEnergiesFile, readStress, readTraj
export getNbStepEnergies, getNbStepStress, getNbStepAtomsFtraj
export writeEnergies, writeStress, writeFtraj
export buildingDataBase

# TODO  Much left to do:
# - Make a cleaning function that deletes the ">>>>>>>>" that occurs
# when doing restarts
# - Check the conversions are ok (CPMD use atomic units )

# Read input file / Extract informations from it
#--------------------------------------------------------------------------------
# Extract timestep of simulation from input file
function getTimestep( file_path::T1 ) where { T1 <: AbstractString }
    # Argument
    # - file_path: path to the input file of CPMD
    # Output
    # - timestep converted to fs
    # OR false if something is wrong with file

    # Check existence of file
    if ! isfile( file_path )
        # If the file doesn't exists, returns false and sends a message
        print( "No input file at ", file_input_path, "\n" )
        return false
    end

    # Read whole input file
    lines = utils.getLines( file_input_path ) # NB: We can probably do better than this...

    # Get number of lines of input file
    nb_lines=size(lines)[1]

    # Initialize timestep to 0
    timestep=0

    # Loop over lines
    for line_nb=1:nb_lines
        # Parse current line with " " deliminator
        keywords = split( lines[line_nb] )

        # Check that the line is not empty
        if size(keywords)[1] > 0
            # Check the TIMESTEP Keyword
            if keywords[1] == "TIMESTEP"
                # The First element of the next line is the timestep value in
                # atomic units, casts string to float
                timestep = parse(Float64, split( lines[line_nb+1] )[1] )
            end
        end
    end

    # Returns timestep and converts it from atomic units to fs
    return conversion.hatime2fs*timestep
end
# Extract the stress tensor stride from CPMD input file
function getStrideStress( file_path::T1 ) where { T1 <: AbstractString }
    # Argument
    # - file_path: path to the input file of CPMD
    # Output
    # - stride_stres: stride of the print of the STRESS file
    # OR false if something went wrong with the file

    # Check that file exists
    if ! isfile( file_path )
        # If not, sends message and returns false
        print("No input file at ",file_input_path,"\n")
        return false
    end

    # Read whole input file
    lines = utils.getLines( file_input_path )

    # Get number of lines of the input file
    nb_lines = size( lines )[1]

    # Initialize stress tensor 0
    stride_stress = 0

    # Loop over lines
    for line_nb=1:nb_lines
        # Parse line with " " deliminator
        keywords=split(lines[line_nb])

        # Check that line is not empty
        if size(keywords)[1] > 0
            # Check that first and second elements of the lines are "STRESS" and "TENSOR"
            if keywords[1] == "STRESS" && keywords[2] == "TENSOR"
                # If so, first element of the next line is the stride
                # Reads it and cast it from String to Float
                stride_stress = parse(Float64, split( lines[line_nb + 1] )[1] )
            end
        end
    end

    # Returns stride of the stress tensor
    return stride_stress
end
# Extract the TRAJ stride from CPMD input file
function getStrideTraj( file_path::T1 ) where { T1 <: AbstractString }
    # Argument
    # - file_path: path to the CPMD input file
    # Output
    # - stride-traj: stride to write the TRAJ file

    # Check that file exists
    if ! isfile( file_path )
        # If not, sends a message and return false
        print("No input file at ",file_input_path,"\n")
        return false
    end

    # Read whole input file
    lines = utils.getLines( file_input_path )

    # Get number of lines of the input file
    nb_lines=size(lines)[1]

    # Initialize stride traj to 0
    stride_traj = 0

    # Loop over lines
    for line_nb=1:nb_lines
        # Parse current line with " " deliminator
        keywords=split(lines[line_nb])

        # Check that line is not empty
        if size(keywords)[1] > 0
            # Check that the first element of the line is "TRAJECTORY"
            if keywords[1] == "TRAJECTORY"
                # If so, the stride of TRAJ is the first element of the next line
                stride_traj = parse(Float64, split( lines[ line_nb + 1 ] )[1] )
            end
        end
    end

    # Return the stride of the trajectory
    return stride_traj
end
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# Write Velocities in the input file (for a soft restart), using IO handler
function writeVelocities( handle_out::T1, velocities::Array{T2,2} ) where { T1 <: IO, T2 <: Real }
    # Arguments
    # - file_out: IO handler of the output file
    # - velocities: velocities of the file (atomic units)
    # Output
    # - Bool, whether writting was successful

    # Write "VELOCITIES" signal
    write( handle_out, string("VELOCITIES\n") )

    # Get number of atoms
    nb_atoms = size(velocities)[1]

    # Write number of atoms
    write( handle_out, string( nb_atoms, " " ) )

    # Loop over atoms
    for atom=1:nb_atoms
        # Write all atomic index
        write( handle_out, string( atom ," ") )
    end

    # Write end of line
    write( handle_out, string("\n") )

    # Loop over atoms
    for atom=1:nb_atoms
        # Loop over dimension
        for i=1:3
            # Write velocities of atom in dimension i
            write( handle_out, string( velocities[atom,i], " " ) )
        end
        # Write end of line for atom
        write( handle_out, string( "\n" ) )
    end

    # Write "END VELOCITIES" signal
    write( handle_out, string("END VELOCITIES\n") )

    # Returns true if we reach this point
    return true
end
# Write Velocities in the input file (for a soft restart), using IO handler
function writeVelocities( handle_out::T1, velocities::Array{T2,2}, atoms_indexes::Vector{T3} ) where { T1 <: IO, T2 <: Real, T3 <: Int }
    # Arguments
    # - handle_out: IO handler for input file
    # - velocities: array (nb_atoms,3) contains the velocities of atoms
    # - atoms_indexes: atoms of the index of the atom to write the velocities in
    # the input file
    # Output
    # - Bool, whether writting was successful

    # Write VELOCITIES signal (start)
    write( handle_out, string("VELOCITIES\n") )

    # Write number of atoms to file
    write( handle_out, string( size(nb_atoms_nb)[1], " " ) )

    # Loop over target atom indexes
    for atom_index in atom_indexes
        # Write atom index
        write( handle_out, string( atom_index, " " ) )
    end

    # Writes end of line
    write( handle_out, string("\n") )

    # Loop over atom indexes
    for atom_index in atoms_indexes
        # Loop over dimensions
        for i=1:3
            # Write velocities of atom_index in dimension i to file
            write( handle_out, string( velocities[atom_index,i], " " ) )
        end
        # Write end of line
        write( handle_out, string( "\n" ) )
    end

    # Write "END VELOCITIES" signal
    write("END VELOCITIES\n")

    # Returns true if we reach this point
    return true
end
# Writes positions to the input file
function writePositions( handle_out::T1, positions::Array{T2,2} ) where { T1 <: IO, T2 <: Real }
    # Argument
    # - handle_out: IO handler of the input file
    # - positions: array (nb_atoms,3) contains positions of the atoms
    # Output
    # - Bool: whether writting was successful

    # Get number of atoms
    nb_atoms=size(positions)[1]

    # Loop over atoms
    for atom=1:nb_atoms
        # Loop over dimensions
        for i=1:3
            # Write atomic position in dimension i
            write( handle_out, string( positions[atom,i], " " ) )
        end

        # Write end of line
        write( handle_out, string( "\n" ) )
    end

    # Returns true if all went ok
    return true
end
#--------------------------------------------------------------------------------

# Compute actual timestep of the output file
#--------------------------------------------------------------------------------
function computeTimestep( file_path::T1 ) where { T1 <: AbstractString }
    # Argument
    # - file_path: path of the input file of CPMD
    # Output
    # - apparent timestep of the output file
    # OR return false if something went wrong

    # Get the timestep of the simulation using input file
    timestep = readInputTimestep( file_path )

    # Check that the timestep reading was ok
    if timestep == false
        return false
    end

    # Get the trajectory stride using input file
    stride_traj = readIntputStrideTraj( file_path )

    # Check that the stride traj reading was ok
    if stride_traj == false
        return false
    end

    # Returns the timestep
    return stride_traj*timestep
end
#--------------------------------------------------------------------------------


# Reading ENERGIES file
#-----------------------------------------
# Contains: Temperature, Potential Energy, Total Energy, MSD and Computing time for each step
# Structure:
# 1 line per step, per column:
# time, temperature, potential energy, total energy, MSD, Computing time
#----------------------------------------------------------------------------
col_time  = 1    # Time elapsed in simulation
col_temp  = 3    # Temperature
col_poten = 4    # Potential Energy
col_entot = 5    # Total Energy
col_msd   = 7    # MSD
col_comp  = 8    # Computing time for SCF
#----------------------------------------------------------------------------
# Get number of step in ENERGIES file
function getNbStepEnergies( file_path::T1 ) where { T1 <: AbstractString }
    # Argument
    # - file_path: path to the ENERGIES file
    # Output
    # - nb_step: number of step in ENERGIES file
    # OR return false if something went wrong

    # Check that the file exists
    if ! isfile( file_path )
        # If not sends a message and returns false
        print("No ENERGIES file at ",file_path,"\n")
        return false
    end

    # Initialize nb_step at 0
    nb_step=0

    # Opens input file
    file_in = open( file_path )

    # Loop over lines, as long as possible
    while ( ! eof(file_in) )
        # readline in empty
        readline( file_in )

        # Increments step counter
        nb_step += 1
    end

    # Number of step
    return nb_step
end
# Reads ENERGIES file
function readEnergies( file_path::T1 ) where { T1 <: AbstractString }
    # Argument
    # - file_path: path to the ENERGIES file
    # Output
    # - temp: Temperature
    # - epot: Potential Energy
    # - etot: Total Energy
    # - msd: Mean Square Displacement
    # - comp: Computational step for SCF
    # OR return false, false, false, false, false

    # Get number of steps in the file
    nb_step = getNbStepEnergies( file_path )

    # If something went wrong reading the number of steps, returns false
    if nb_step == false
        return false, false, false, false, false
    end

    # Initialize data vector
    temp = zeros( nb_step ) # Temperature
    epot = zeros( nb_step ) # Potential Energy
    etot = zeros( nb_step ) # Total Energy
    msd  = zeros( nb_step ) # MSD
    comp = zeros( nb_step ) # Computational time for SCF

    # Opens input file (ENERGIES)
    handle_in = open(file_path )

    # Loop over steps
    for step=1:nb_step
        # Read line and parse with " " deliminator
        line = split( readline(file_in) )

        # Parse data and attribute data to each vector
        temp[step] = parse(Float64, line[ col_temp  ] ) # Temp
        epot[step] = parse(Float64, line[ col_poten ] ) # Potential Energy
        etot[step] = parse(Float64, line[ col_entot ] ) # Total energy
        msd[step]  = parse(Float64, line[ col_msd   ] ) # MSD
        comp[step] = parse(Float64, line[ col_comp  ] ) # Computational time
    end

    # Close input file
    close(file_in)

    # Returns the vectors with data
    return  temp, epot, etot, msd, comp
end
# Reads ENERGIES file with a given stride
function readEnergies( file_path::T1, stride_::T2 ) where { T1 <: AbstractString, T2 <: Int  }

    # Check file
    #----------------------------------------
    nb_step_origin =  getNbStepEnergies( file_path )
    if nb_step_origin == false
        return false, false, false, false, false
    end
    #----------------------------------------

    # Array Init
    #----------------------------------------
    nb_step = utils.nbStepStriding( nb_step_origin, stride_ )
    temp = zeros( nb_step )
    epot = zeros( nb_step )
    etot = zeros( nb_step )
    msd  = zeros( nb_step )
    comp = zeros( nb_step )
    #----------------------------------------

    # Getting data from lines
    #----------------------------------------------
    file_in = open(file_path)
    count_=1
    for step=1:nb_step_origin
        line=split( readline(file_in) )
        if (step-1) % stride_ == 0
            temp[count_] = parse( Float64, line[col_temp] )
            epot[count_] = parse( Float64, line[col_poten] )
            etot[count_]  = parse( Float64, line[col_entot] )
            msd[count_]   = parse( Float64, line[col_msd] )
            comp[count_]  = parse( Float64, line[col_comp] )
            count_ += 1
        end
    end
    close(file_in)
    #----------------------------------------------

    return  temp, epot, etot, msd, comp
end
function readEnergies( file_path::T1, stride_::T2, nb_ignored::T3 ) where { T1 <: AbstractString, T2 <: Int, T3 <: Int  }

    # Check file
    #----------------------------------------
    nb_step_origin =  getNbStepEnergies( file_path )
    if nb_step_origin == false
        return false, false, false, false, false
    end
    #----------------------------------------

    # Array Init
    #----------------------------------------
    nb_step = utils.nbStepStriding( nb_step_origin-nb_ignored, stride_ )
    temp = zeros(nb_step)
    epot = zeros(nb_step)
    etot = zeros(nb_step)
    msd  = zeros(nb_step)
    comp = zeros(nb_step)
    #----------------------------------------

    # Getting data from lines
    #----------------------------------------------
    file_in = open(file_path)
    count_=1
    for step=1:nb_ignored
        temp=readline(file_in)
    end
    for step=1:nb_step_origin-nb_ignored
        line=split( readline(file_in) )
        if (step-1) % stride_ == 0 && step
            temp[count_] = parse( Float64, line[col_temp] )
            epot[count_] = parse( Float64, line[col_poten] )
            etot[count_]  = parse( Float64, line[col_entot] )
            msd[count_]   = parse( Float64, line[col_msd] )
            comp[count_]  = parse( Float64, line[col_comp] )
            count_ += 1
        end
    end
    close(file_in)
    #----------------------------------------------

    return  temp, epot, etot, msd, comp
end
function readEnergies( file_path::T1, stride_::T2, nb_ignored::T3, nb_max::T4 ) where { T1 <: AbstractString, T2 <: Int, T3 <: Int, T4 <: Int  }

    # Check file
    nb_step_origin =  getNbStepEnergies( file_path )
    if nb_step_origin == false
        return false, false, false, false, false
    end

    # Array Init
    #---------------------------------------
    nb_step = utils.nbStepStriding( nb_step_origin-nb_ignored, stride_ )
    if nb_max > nb_step
        print("nb_max is too large, maximum value is ",nb_step,"\n")
    end
    if nb_max <= 0
        print("nb_max must be positive!\n")
    end
    temp=zeros(nb_max)
    epot=zeros(nb_max)
    etot=zeros(nb_max)
    msd=zeros(nb_max)
    comp=zeros(nb_max)
    #----------------------------------------

    # Getting data from lines
    #----------------------------------------------
    file_in = open(file_path)
    for step=1:nb_ignored
        temp=readline(file_in)
    end
    count_=1
    for step=1:nb_step_origin - nb_ignored
        line=split( readline(file_in) )
        if (step-1) % stride_ == 0
            temp[count_] = parse( Float64, line[col_temp] )
            epot[count_] = parse( Float64, line[col_poten] )
            etot[count_]  = parse( Float64, line[col_entot] )
            msd[count_]   = parse( Float64, line[col_msd] )
            comp[count_]  = parse( Float64, line[col_comp] )
            if count_ >= nb_max
                break
            end
            count_ += 1
        end
    end
    close(file_in)
    #----------------------------------------------

    return  temp, epot, etot, msd, comp
end
function writeEnergies( file_path::T1, temp::Vector{T2}, epot::Vector{T3}, etot::Vector{T4}, msd::Vector{T5}, comp::Vector{T6} ) where { T1 <: AbstractString, T2 <: Real, T3 <: Real, T4 <: Real, T5 <: Real, T6 <: Real }
    nb_step=size(temp)[1]
    file_output = open( file_path, "w" )
    for step=1:nb_step
        write( file_output, string( step, " " ) )
        write( file_output, string( 0, " " ) )
        write( file_output, string( temp[step], " " ) )
        write( file_output, string( epot[step], " " ) )
        write( file_output, string( etot[step], " " ) )
        write( file_output, string( 0, " " ) )
        write( file_output, string( msd[step], " " ) )
        write( file_output, string( comp[step], " " ) )
        write( file_output, string( "\n" ) )
    end
    close(file_output)
    return true
end
#----------------------------------------------------------------------------


#------------------------
# Reading STRESS file
#------------------------------------------------------------------------------#
# Contains: Stress tensor for each step (with a possible stride)
# Structure:
# 4 lines per step
# line 1: Comment (indicates step number)
# line 2-4: stress tensor in matrix form
# Sxx Sxy Sxz
# Syx Syy Syz
# Szx Szy Szz
#------------------------------------------------------------------------------#
stress_block_size=4
stress_dim=3
function getNbStepStress( file_path::T1 ) where { T1 <: AbstractString }
    # Check if file exists
    if ! isfile( file_path )
        print("STRESS file does not exists at ",file_path," !\n")
        return false
    end
    # Counting the Number of lines
    nb_line=0
    file_in = open( file_path )
    while ( ! eof(file_in) )
        temp=readline(file_in)
        nb_line += 1
    end
    close(file_in)
    # If the number of lines is not nb_line*block_size, the file is likely corrupted
    if nb_line % stress_block_size != 0
        print("STRESS file likely corrupted at ",file_path," !\n")
        return false
    end
    # Returns number of blocks
    return Int(nb_line/stress_block_size)
end
function readStress( file_path::T1 ) where { T1 <: AbstractString }

    # Checking file exists
    nb_step = getNbStepStress( file_path )
    if nb_step == false
        return false
    end

    # Init data files
    #------------------------------------------
    stress = zeros( Real, nb_step, stress_dim, stress_dim )
    #------------------------------------------

    #--------------------------------------------------------
    file_in = open( file_path )
    for step=1:nb_step
        temp=split( readline( file_in ) ) # Comment line
        if temp[1] != "TOTAL"
            # Corruption in file
            print("STRESS file is likely corrupted!\n")
            print("line: ",step*stress_block_size,"\n")
            return false
        end
        for i=1:stress_dim
            keywords=split( readline( file_in ) )
            if keywords[1] == "TOTAL"
                print("STRESS file is likely corrupted!\n")
                print("line: ",step*stress_block_size+i,"\n")
                return false
            else
                for j=1:stress_dim
                    stress[step,i,j] = parse(Float64,keywords[j])
                end
            end
        end
    end
    close(file_in)
    #--------------------------------------------------------

    return stress
end
function readStress( file_path::T1, stride_::T2 ) where { T1 <: AbstractString, T2 <: Int }

    # Checking file exists
    nb_step_origin = getNbStepStress( file_path )
    if nb_step_origin == false
        return false
    end

    # Init data files
    #------------------------------------------
    nb_step = utils.nbStepStriding( nb_step_origin, stride_ )
    stress  = zeros( Real, nb_step, stress_dim, stress_dim )
    #------------------------------------------

    #--------------------------------------------------------
    count_=1
    file_in=open(file_path)
    for step=1:nb_step_origin
        if (step-1) % stride_ == 0
            temp=split( readline( file_in ) ) # Comment line
            if temp[1] != "TOTAL"
                # Corruption in file
                print("STRESS file is likely corrupted!\n")
                print("line: ",step*stress_block_size,"\n")
                return false
            end
            for i=1:stress_dim
                keywords=split( readline( file_in ) )
                if keywords[1] == "TOTAL"
                    print("STRESS file is likely corrupted!\n")
                    print("line: ",step*stress_block_size+i,"\n")
                    return false
                else
                    for j=1:stress_dim
                        stress[count_,i,j] = parse(Float64,keywords[j])
                    end
                end
            end
            count_ += 1
        else
            for skip_line=1:stress_block_size
                temp = readline( file_in )
            end
        end
    end
    close(file_in)
    #--------------------------------------------------------

    return stress
end
function readStress( file_path::T1, stride_::T2, nb_ignored::T3 ) where { T1 <: AbstractString, T2 <: Int, T3 <: Int }

    # Checking file exists
    nb_step_origin = getNbStepStress( file_path )
    if nb_step_origin == false
        return false
    end

    # Init data files
    #------------------------------------------
    nb_step = utils.nbStepStriding( nb_step_origin-nb_ignored, stride_ )
    stress = zeros( Real, nb_step, stress_dim, stress_dim )
    #------------------------------------------

    #--------------------------------------------------------
    file_in=open(file_path)
    for step=1:nb_ignored
        for i=1:stress_block_size
            temp = readline(file_in )
        end
    end
    count_=1
    for step=1:(nb_step_origin-nb_ignored)
        if (step-1) % stride_ == 0
            temp=split( readline( file_in ) ) # Comment line
            if temp[1] != "TOTAL"
                # Corruption in file
                print("STRESS file is likely corrupted!\n")
                print("line: ",step*stress_block_size,"\n")
                return false
            end
            for i=1:stress_dim
                keywords=split( readline( file_in ) )
                if keywords[1] == "TOTAL"
                    print("STRESS file is likely corrupted!\n")
                    print("line: ",step*stress_block_size+i,"\n")
                    return false
                else
                    for j=1:stress_dim
                        stress[count_,i,j] = parse(Float64,keywords[j])
                    end
                end
            end
            count_+=1
        else
            for skip_line=1:stress_block_size
                temp = readline( file_in )
            end
        end
    end
    close(file_in)
    #--------------------------------------------------------

    return stress
end
function readStress( file_path::T1, stride_::T2, nb_ignored::T3, nb_max::T4 ) where { T1 <: AbstractString, T2 <: Int, T3 <: Int, T4 <: Int }

    nb_step_origin=getNbStepStress( file_path )
    if nb_step_origin == false
        return false
    end

    # Init data files
    #------------------------------------------
    nb_step = utils.nbStepStriding( nb_step_origin-nb_ignored, stride_ )
    if nb_max > nb_step
        print("nb_max is too large, maximum value is ",nb_step,"\n")
    end
    if nb_max <= 0
        print("nb_max must be positive!\n")
    end
    stress=zeros( Real, nb_max, stress_dim, stress_dim )
    #------------------------------------------

    #--------------------------------------------------------
    file_in=open(file_path)
    for step=1:nb_ignored
        for i=1:stress_block_size
            temp = readline(file_in )
        end
    end
    count_=1
    for step=1:(nb_step_origin-nb_ignored)
        if (step-1) % stride_ == 0
            temp=split( readline( file_in ) ) # Comment line
            if temp[1] != "TOTAL"
                # Corruption in file
                print("STRESS file is likely corrupted!\n")
                print("line: ",step*stress_block_size,"\n")
                return false
            end
            for i=1:stress_dim
                keywords=split( readline( file_in ) )
                if keywords[1] == "TOTAL"
                    print("STRESS file is likely corrupted!\n")
                    print("line: ",step*stress_block_size+i,"\n")
                    return false
                else
                    for j=1:stress_dim
                        stress[count_,i,j] = parse(Float64,keywords[j])
                    end
                end
            end
            if count_ >= nb_max
                break
            end
            count_+=1
        else
            for skip_line=1:stress_block_size
                temp = readline( file_in )
            end
        end
    end
    close(file_in)
    #--------------------------------------------------------

    return stress
end
function writeStress( file_path::T1, stress_tensor::Array{T2,3} ) where { T1 <: AbstractString, T2 <: Real }
    nb_step=size(stress_tensor)[1]
    file_out = open( file_path, "w" )
    for step=1:nb_step
        write(file_out,string("TOTAL STRESS TENSOR (kB): STEP: ",step,"\n"))
        for i=1:3
            for j=1:3
                write(file_out,string(stress_tensor[step,i,j]," "))
            end
            write(file_out,string("\n"))
        end
    end
    close(file_out)
    return true
end
#------------------------------------------------------------------------------#


# Read FTRAJECTORY file
#-------------------------------------------------------------------------------
# Contains: positions, velocity and forces in atomic units for each step
# Structure:
# 1 line per atom per step
# Per line:
# atom_number time x y z vx vy vz fx fy fz
#-------------------------------------
# Time in timestep (accounts for stride_traj)
# Positions in Bohr
# Velocities in Bohr/tHart
# Forces in Ha/Bohr
#--------------------------------------
# When restarting run, a misc line that needs to be ignored
# <<<<<<  NEW DATA  >>>>>>
#-------------------------------------
col_start_velocity=4
col_start_position=1
col_start_force=7
function getNbStepAtomsFtraj( file_path::T1 ) where { T1 <: AbstractString }

    if ! isfile( file_path )
        print("File FTRAJECTORY does not exists at ",file_path,"\n")
        return false, false
    end

    file_in = open( file_path )
    nb_line = 0
    nb_atoms=0
    while ! eof( file_in )
        keywords=split(readline( file_in ))
        if keywords[1] != "<<<<<<"
            if keywords[1] == "1"
                nb_atoms += 1
            end
            nb_line += 1
        end
    end
    close( file_in )

    if nb_line % nb_atoms != 0
        print("Potential Corruption in file FTRAJECTORY at ",file_path,"\n")
        print("Number of line: ", nb_line, "\n" )
        print("Number of atoms: ", nb_atoms, "\n" )
        print("Modulo: ", nb_line%nb_atoms, "\n" )
        return false, false
    end

    return Int(nb_line/nb_atoms), nb_atoms
end
function readFtraj( file_path::T1 ) where { T1 <: AbstractString }

    # Getting number of line of file and checking existence of file
    #-----------------------------------------------------
    nb_step, nb_atoms = getNbStepAtomsFtraj( file_path )
    if nb_step == false
        return false, false, false
    end
    #-----------------------------------------------------

    # Init arrays
    #-----------------------------------------------------
    positions  = zeros( Real, nb_step, nb_atoms, 3 )
    velocities = zeros( Real, nb_step, nb_atoms, 3 )
    forces     = zeros( Real, nb_step, nb_atoms, 3 )
    #-----------------------------------------------------

    # Reading
    #-----------------------------------------------------
    file_in = open( file_path )
    for step=1:nb_step
        for atom=1:nb_atoms
            keywords=split( readline( file_in ) )
            if keywords[1] != "<<<<<<"
                for i=1:3
                    positions[step,atom,i]  = parse(Float64, keywords[ i + col_start_position ] )
                    velocities[step,atom,i] = parse(Float64, keywords[ i + col_start_velocity ] )
                    forces[step,atom,i]     = parse(Float64, keywords[ i + col_start_force ] )
                end
            end
        end
    end
    close( file_in )
    #-----------------------------------------------------

    return positions, velocities, forces
end
function readFtraj( file_path::T1, stride_::T2 ) where { T1 <: AbstractString, T2 <: Int }

    # Getting number of line of file
    #---------------------------------------------------------------------------
    nb_step_origin, nb_atoms = getNbStepAtomsFtraj( file_path )
    if nb_step_origin == false
        print("File FTRAJ does not exists!\n")
        return false, false, false
    end
    #---------------------------------------------------------------------------

    #---------------------------------------------------------------------------
    nb_step = utils.nbStepStriding( nb_step_origin, stride_ )
    positions  = zeros( Real, nb_step, nb_atoms, 3 )
    velocities = zeros( Real, nb_step, nb_atoms, 3 )
    forces     = zeros( Real, nb_step, nb_atoms, 3 )
    #---------------------------------------------------------------------------

    #---------------------------------------------------------------------------
    file_in = open( file_path )
    count_step=1
    for step=1:nb_step_origin
        if (step-1) % stride_ == 0
            for atom=1:nb_atoms
                keywords=split( readline( file_in ) )
                if keywords[1] == "<<<<<<"
                    keywords=split( readline( file_in ) )
                end
                for i=1:3
                    positions[count_step,atom,i]  = parse(Float64, keywords[ i + col_start_position ] )
                    velocities[count_step,atom,i] = parse(Float64, keywords[ i + col_start_velocity ] )
                    forces[count_step,atom,i]     = parse(Float64, keywords[ i + col_start_force ] )
                end
            end
            count_step += 1
        end
    end
    close( file_in )
    #---------------------------------------------------------------------------

    return positions, velocities, forces
end
function readFtraj( file_path::T1, stride_::T2, nb_ignored::T3 ) where { T1 <: AbstractString, T2 <: Int, T3 <: Int }

    # Getting number of line of file
    #---------------------------------------------------------------------------
    nb_step_origin, nb_atoms = getNbStepAtomsFtraj( file_path )
    if nb_step_origin == false
        print("File FTRAJ does not exists!\n")
        return false, false, false
    end
    #---------------------------------------------------------------------------

    #---------------------------------------------------------------------------
    nb_step = utils.nbStepStriding( nb_step_origin-nb_ignored, stride_ )
    positions  = zeros( Real, nb_step, nb_atoms, 3 )
    velocities = zeros( Real, nb_step, nb_atoms, 3 )
    forces     = zeros( Real, nb_step, nb_atoms, 3 )
    #---------------------------------------------------------------------------

    #---------------------------------------------------------------------------
    file_in = open( file_path )
    # Ignoring the first nb_ignore steps
    for step=1:nb_ignored
        for atom=1:nb_atoms
            temp=readline(file_in)
        end
    end
    count_step=1
    for step=1:nb_step_origin-nb_ignored
        if (step-1) % stride_ == 0
            for atom=1:nb_atoms
                keywords=split( readline( file_in ) )
                if keywords[1] == "<<<<<<"
                    keywords=split( readline( file_in ) )
                end
                for i=1:3
                    positions[count_step,atom,i]  = parse(Float64, keywords[ i + col_start_position ] )
                    velocities[count_step,atom,i] = parse(Float64, keywords[ i + col_start_velocity ] )
                    forces[count_step,atom,i]     = parse(Float64, keywords[ i + col_start_force ] )
                end
            end
            count_step += 1
        end
    end
    close( file_in )
    #---------------------------------------------------------------------------

    return positions, velocities, forces
end
function readFtraj( file_path::T1, stride_::T2, nb_ignored::T3, nb_max::T4 ) where { T1 <: AbstractString, T2 <: Int, T3 <: Int, T4 <: Int }

    # Getting number of line of file
    #---------------------------------------------------------------------------
    nb_step_origin, nb_atoms = getNbStepAtomsFtraj( file_path )
    if nb_step_origin == false
        print("File FTRAJ does not exists!\n")
        return false, false, false
    end
    #---------------------------------------------------------------------------

    #---------------------------------------------------------------------------
    nb_step = utils.nbStepStriding( nb_step_origin-nb_ignored, stride_ )
    if nb_max > nb_step
        print("nb_max is too large, maximum value is ",nb_step,"\n")
    end
    if nb_max <= 0
        print("nb_max must be positive!\n")
    end
    positions  = zeros( nb_max, nb_atoms, 3 )
    velocities = zeros( nb_max, nb_atoms, 3 )
    forces     = zeros( nb_max, nb_atoms, 3 )
    #---------------------------------------------------------------------------

    #---------------------------------------------------------------------------
    file_in = open( file_path )
    # Ignoring the first nb_ignore steps
    for step=1:nb_ignored
        for atom=1:nb_atoms
            temp=readline(file_in)
        end
    end
    count_step=1
    for step=1:nb_step_origin-nb_ignored
        if (step-1) % stride_ == 0
            for atom=1:nb_atoms
                keywords=split( readline( file_in ) )
                if keywords[1] == "<<<<<<"
                    keywords=split( readline( file_in ) )
                end
                for i=1:3
                    positions[count_step,atom,i]  = parse(Float64, keywords[ i + col_start_position ] )
                    velocities[count_step,atom,i] = parse(Float64, keywords[ i + col_start_velocity ] )
                    forces[count_step,atom,i]     = parse(Float64, keywords[ i + col_start_force ] )
                end
            end
            if count_step >= nb_max
                break
            end
            count_step += 1
        end
    end
    close( file_in )
    #---------------------------------------------------------------------------

    return positions, velocities, forces
end
function writeFtraj( file_path::T1, positions::Array{T2,3}, velocities::Array{T3,3}, forces::Array{T4,3} ) where { T1 <: AbstractString, T2 <: Real, T3 <: Real, T4 <: Real }
    nb_step = size(positions)[1]
    nb_atoms = size(positions)[2]
    file_out = open( file_path, "w" )
    for step=1:nb_step
        for atom=1:nb_atoms
            write( file_out, string( step, " " ) )
            for i=1:3
                write( file_out, string( positions[step,atom,i], " " ) )
            end
            for i=1:3
                write( file_out, string( velocities[step,atom,i], " " ) )
            end
            for i=1:3
                write( file_out, string( forces[step,atom,i], " " ) )
            end
            write(file_out,string("\n"))
        end
    end
    close( file_out )
    return true
end
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
function buildingDataBase( folder_target::T1, file_stress::T2, file_pressure::T3, file_traj::T4, file_ftraj::T5, file_energy::T6, timestep_target::T7 ) where { T1 <: AbstractString, T2 <: AbstractString, T3 <: AbstractString, T4 <: AbstractString, T5 <: AbstractString, T6 <: AbstractString, T7 <: Real }

    # Reading input file
    #---------------------------------------------------------------------------
    file_input=string(folder_target,"input")
    # Getting the stride for the STRESS file
    stride_stress = cpmd.readIntputStrideStress( file_input )
    # Getting the stride for the TRAJEC and FTRAJECTORY files
    stride_traj   = cpmd.readIntputStrideTraj( file_input )
    # Getting the timestep of the simulation
    timestep_sim  =  cpmd.readInputTimestep( file_input )
    #---------------------------------------------------------------------------

    # Computing
    #---------------------------------------------------------------------------
    n_stress = round( Int, timestep_target/( timestep_sim*stride_stress ) )
    n_traj   = round( Int, timestep_target/( timestep_sim*stride_traj ) )
    n_energy = round( Int, timestep_target/( timestep_sim ) )
    n_ftraj  = round( Int, timestep_target/( timestep_sim*stride_traj ) )
    #---------------------------------------------------------------------------

    # Target files
    #---------------------------------------------------------------------------
    file_energy_in = string(folder_target,"ENERGIES")
    file_trajec_in = string(folder_target,"TRAJEC.xyz")
    file_stress_in = string(folder_target,"STRESS")
    file_ftrajectory_in = string(folder_target,"FTRAJECTORY")
    #---------------------------------------------------------------------------

    # Determining nb of steps
    #---------------------------------------------------------------------------
    nb_step_stress = getNbStepStress( file_stress_in )
    if nb_step_stress == false
        print("Interrupting database construction.\n")
        return false
    end
    nb_step_stress = utils.nbStepStriding( nb_step_stress, n_stress )
    #-------------------------------------------
    nb_step_ftraj, nb_atoms_ftraj = getNbStepAtomsFtraj( file_ftrajectory_in )
    if nb_step_ftraj == false
        print("Interrupting database construction.\n")
        return false
    end
    nb_step_ftraj = utils.nbStepStriding( nb_step_ftraj, n_ftraj )
    #-------------------------------------------
    nb_step_energy = getNbStepEnergies( file_energy_in )
    if nb_step_energy == false
        print("Interrupting database construction.\n")
        return false
    end
    nb_step_energy = utils.nbStepStriding( nb_step_energy, n_energy )
    #-------------------------------------------
    nb_step_traj = filexyz.getNbSteps( file_trajec_in )
    if nb_step_traj == false
        print("Interrupting database construction.\n")
        return false
    end
    nb_step_traj = utils.nbStepStriding( nb_step_traj, n_traj )
    #--------------------------------------------------------------------------

    # Checking coherence
    #--------------------------------------------------------------------------
    if nb_step_traj != nb_step_stress ||  nb_step_traj != nb_step_energy || nb_step_traj != nb_step_ftraj
        print("Some inconsistencies in ",folder_target,"\n")
        print("traj_step: ",nb_step_traj,"\n")
        print("ftraj_step: ",nb_step_ftraj,"\n")
        print("energy_step: ",nb_step_energy,"\n")
        print("stress_step: ",nb_step_stress,"\n")
        print("Interrupting database construction.\n")
        return false
    end
    target_length = min( nb_step_traj, nb_step_ftraj, nb_step_energy, nb_step_stress )
    #--------------------------------------------------------------------------

    # Writing
    #--------------------------------------------------------------------------
    nb_ignored=0
    stress_tensor = readStress( file_stress_in, n_stress, nb_ignored, target_length )
    writeStress( file_stress, stress_tensor )
    utils.writeData( file_pressure, press_stress.computePressure(stress_tensor) )
    stress_tensor=[] # Clearing memory
    filexyz.writeXYZ( file_traj, filexyz.readFileAtomList( file_trajec_in, n_traj, nb_ignored, target_length ) )
    positions,velocities,forces=readFtraj( file_ftrajectory_in, n_ftraj, nb_ignored, target_length )
    writeFtraj( file_ftraj, positions, velocities, forces )
    temp, epot, etot, msd, comp = readEnergies( file_energy_in, n_energy, nb_ignored, target_length )
    writeEnergies( file_energy, temp, epot, etot, msd, comp )
    #--------------------------------------------------------------------------

    return true
end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
function writeRelaunchPositions( folder_target::T1, atoms::T2, file_out_positions_suffix::T3 ) where { T1 <: AbstractString, T2 <: atom_mod.AtomList, T3 <: AbstractString }

    #--------------------------------------------------------
    species = atom_mod.getSpecies( atoms )
    nb_element_species = atom_mod.getNbElementSpecies( atoms )
    start_species = atom_mod.getStartSpecies( atoms )
    nb_species=size(species)[1]
    #--------------------------------------------------------

    #--------------------------------------------------------
    for i_spec=1:nb_species
        file_out=open( string( folder_in_target, species[i_spec], "_", file_out_positions_suffix ), "w" )
        for atom = start_species[i_spec]:start_species[i_spec]+nb_element_species[i_spec]
            for i=1:3
                write( file_out, string( atoms.positions[atom,i]), " " )
            end
            write(file_out,"\n")
        end
        close( file_out )
    end
    #--------------------------------------------------------

    return true
end
function writeRelaunchVelocities( folder_target::T1, velocities::Array{T2,2}, file_out_velocities::T3 ) where { T1 <: AbstractString, T2 <: Real, T3 <: AbstractString }

    #----------------------------------------
    nb_atoms=size(velocities)[1]
    file_out = open( string( folder_target, file_out_velocities ), "w" )
    for atom = 1:nb_atoms
        for i=1:3
            write( file_out, string( velocities[atom,i], " " ) )
        end
        write("\n")
    end
    close( file_out )
    #----------------------------------------

    return true
end
function writeRelaunchVelocitiesTraj( folder_target::T1, traj::Vector{T2}, file_out_velocities::T3 ) where { T1 <: AbstractString, T2 <: atom_mod.AtomList, T3 <: AbstractString }

    #----------------------------------------
    nb_step  = filexyz.getNbStep( traj )
    nb_atoms = filexyz.getNbAtoms( traj[nb_step] )
    #----------------------------------------

    #---------------------------------------
    timestep = readInputTimestep( string( folder_target, "input" ) )
    if timestep == false
        return false
    end
    stride_traj = readIntputStrideTraj( string( folder_target, "input" ) )
    if stride_traj == false
        return false
    end
    dt = stride_traj*timestep*conversion.fs2hatime
    #---------------------------------------

    #----------------------------------------
    file_out = open( string( folder_target, file_out_velocities ), "w" )
    for atom = 1:nb_atoms
        for i=1:3
            dx=(traj[nb_step].positions[atom,i]-traj[nb_step-1].positions[atom,i])*conversion.ang2Bohr
            write( file_out, string( dx/dt , " " ) )
        end
        write("\n")
    end
    close( file_out )
    #----------------------------------------

    return true
end
function relaunchRunTrajec( folder_in_target::T1, file_out_positions_suffix::T2, file_out_velocities::T3 ) where { T1 <: AbstractString, T2 <: AbstractString, T3 <: AbstractString }

    #------------------------------------------------
    traj = filexyz.readFileAtomList( string( folder_in_target, "TRAJEC.xyz" ) )
    if traj == false
        return false
    end
    nb_step_traj = size(traj)[1]
    #------------------------------------------------

    #------------------------------------------------
    if writeRelaunchPositions( folder_in_target, traj[1], file_out_positions_suffix ) && writeRelaunchVelocitiesTraj( folder_in_target, traj, file_out_velocities )
        return true
    else
        return false
    end
    #-------------------------------------------------
end
function relaunchRunFtraj( folder_in_target::T1, file_out_positions_suffix::T2, file_out_velocities::T3 ) where { T1 <: AbstractString, T2 <: AbstractString, T3 <: AbstractString }

    #------------------------------------------------
    positions, velocities, forces = readFTraj( string( folder_in_target, "FTRAJECTORY" ) )
    if positions == false
        return false
    end
    #------------------------------------------------

    #------------------------------------------------
    nb_step_ftraj = size( positions )[1]
    positions = positions[nb_step_ftraj,:,:]
    velocities = velocities[nb_step_ftraj,:,:]
    forces=[]
    #------------------------------------------------

    #------------------------------------------------
    traj = filexyz.readFileAtomList( string( folder_in_target, "TRAJEC.xyz" ) )
    nb_step_traj = size(traj)[1]
    #------------------------------------------------

    #------------------------------------------------
    if nb_step_traj != nb_step_ftraj
        print("TRAJEC.xyz and FTRAJECTORY DO NOT HAVE THE SAME SIZE!\n")
        print("TRAJEC.xyz at: ",folder_in_target,"TRAJEC.xyz n_step: ",nb_step_traj,"\n")
        print("FTRAJECTORY at: ",folder_in_target,"FTRAJECTORY n_step: ",nb_step_ftraj,"\n")
        return false
    end
    #------------------------------------------------

    #------------------------------------------------
    if writeRelaunchPositions( folder_in_target, traj[nb_step_traj], file_out_positions_suffix ) && writeRelaunchVelocities( folder_in_target, velocities, file_out_velocities )
        return true
    end
end
function relaunchRunTrajec( folder_in_target::T1, file_out_path::T2 ) where { T1 <: AbstractString, T2 <: AbstractString }

    #------------------------------------------------------------
    path_input_file = string(folder_in_target,"input")
    if ! isfile(path_input_file)
        print("No input file found at ",path_input_file,"!\n")
        return false
    end
    file_in=open( path_input_file )
    #------------------------------------------------------------

    # Getting positions and velocities
    #------------------------------------------------------------
    dt = computeTimestep( path_input_file )
    if dt == false
        return false
    end
    traj = filexyz.readFileAtomList( string( folder_in_target, "TRAJEC.xyz" ) )
    if traj == false
        return false
    end
    nb_step=size(traj)[1]
    velocities = atom_mod.computeVelocities( traj, nb_step, dt )*conversion.ang2Bohr/conversion.fs2hatime
    traj=traj[nb_step]
    #------------------------------------------------------------

    # Copy parameters of input
    #---------------------------------------
    file_out = open( file_out_path, "w")
    copyInputParams( file_in, file_out )
    #----------------------------------------
    # Copying &ATOMS line
    write(file_out,"&ATOMS\n")
    # Looping over atoms species
    while true
        keywords  = utils.getLineElements( file_in )
        # Skip blank lines
        while size(keywords)[1] == 0
            keywords  = utils.getLineElements( file_in )
        end
        if keywords[1] == "&END" || keywords[1] == "VELOCITIES"
            break
        end
        # Copy PP line
        utils.copyLine2file( keywords, file_out )
        # Copy basis PP line
        utils.copyLine2file( utils.getLineElements( file_in ), file_out )
        # Getting Nb of atoms of species
        keywords =  utils.getLineElements( file_in )
        nb_atoms = parse( Int, keywords[1] )
        # Writing nb atoms line to file
        utils.copyLine2file( keywords, file_out )
        # Writting actual positions for specie
        writePositions( file_out, traj.positions )
        # Ignore the atoms positions
        utils.skipLines( file_in, nb_atoms)
    end
    write( file_out, string("\n") )
    writeVelocities( file_out, velocities )
    write( file_out, string("&END") )
    close( file_in  )
    close( file_out )
end
function relaunchRunFtraj( folder_in_target::T1, file_out_path::T2 ) where { T1 <: AbstractString, T2 <: AbstractString }

    #------------------------------------------------------------
    path_input_file = string(folder_in_target,"input")
    if ! isfile(path_input_file)
        print("No input file found at ",path_input_file,"!\n")
        return false
    end
    file_in=open( path_input_file )
    #------------------------------------------------------------

    # Getting positions and velocities
    #------------------------------------------------------------
    positions, velocities, forces = readFtraj( string( folder_in_target, "FTRAJECTORY_db" ) )
    if positions == false
        return false
    end
    forces = []
    nb_step = size( positions )[1]
    #------------------------------------------------------------

    # Copy parameters of input
    #---------------------------------------
    file_out = open( file_out_path, "w")
    copyInputParams( file_in, file_out )
    #----------------------------------------
    # Copying &ATOMS line
    atom_done=1
    write(file_out,"&ATOMS\n")
    # Looping over atoms species
    while true
        keywords  = utils.getLineElements( file_in )
        # Skip blank lines
        while size(keywords)[1] == 0
            keywords  = utils.getLineElements( file_in )
        end
        if keywords[1] == "&END" || keywords[1] == "VELOCITIES"
            break
        end
        # Copy PP line
        utils.copyLine2file( keywords, file_out )
        # Copy basis PP line
        utils.copyLine2file( utils.getLineElements( file_in ), file_out )
        # Getting Nb of atoms of species
        keywords =  utils.getLineElements( file_in )
        nb_specie = parse( Int, keywords[1] )
        # Writing nb atoms line to file
        utils.copyLine2file( keywords, file_out )
        # Writting actual positions for specie
        writePositions( file_out, positions[nb_step,atom_done:atom_done+nb_specie-1,:]*conversion.bohr2Ang )
        atom_done += nb_specie
        # Ignore the atoms positions
        utils.skipLines( file_in, nb_specie )
    end
    write( file_out, string("\n") )
    writeVelocities( file_out, velocities[nb_step,:,:] )
    write( file_out, string("&END") )
    close( file_in  )
    close( file_out )
end
#-------------------------------------------------------------------------------

function cpmdArcheology( folder_target::T1 , output_suffix::T2 ) where { T1 <: AbstractString, T2 <: AbstractString }

    #---------------------------------------------------------
    if ! isdir( folder_target )
        print("No folder at ",folder_target," !\n")
        return false
    end
    if ! isfile( string( folder_target,"TRAJEC.xyz"))
        print("No TRAJEC.xyz file at ",folder_target," , this function will not work, try another one.\n")
        return false
    end
    if ! isfile( string( folder_target,"ENERGIES") )
        print("No ENERGIES file in ",folder_target," , this function will not work, try another one.\n")
        return false
    end
    if ! isfile( string( folder_target, "FTRAJECTORY" ) )
        print("No FTRAJECTORY file in ", folder_target, " this function will not work, try another one.\n")
        return false
    end
    if ! isifile( string( folder_target, "STRESS" ) )
        print("No STRESS file in ", folder_target, " this function will not work, try another one.\n")
        return false
    end
    if ! isfile( string( folder_target, "input") )
        print("No input file at ",folder_target," \n")
        return false
    end
    #---------------------------------------------------------

    #---------------------------------------------------------
    file_in_traj    = open( string( folder_target, "TRAJEC.xyz" ) )
    file_in_stress  = open( string( folder_target, "STRESS") )
    file_in_ftraj   = open( string( folder_target, "FTRAJECTORY" ) )
    file_in_energy  = open( string( folder_target, "ENERGIES" ) )
    #---------------------------------------------------------

    #---------------------------------------------------------
    file_out_traj   = open( string( folder_target, "TRAJEC", output_suffix, ".xyz" ), "w" )
    file_out_ftraj  = open( string( folder_target, "FTRAJECTORY", output_suffix ), "w" )
    file_out_stress = open( string( folder_target, "STRESS", output_suffix ), "w" )
    file_out_energy = open( string( folder_target, "ENERGIES", output_suffix ), "w" )
    #---------------------------------------------------------

    # Reading steps
    step_number=1

    return true
end

end
