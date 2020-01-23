module cpmd

using conversion
using utils
using filexyz
using press_stress
using atom_mod

export readInputTimestep, readIntputStrideStress, readIntputStrideTraj
export readEnergiesFile, readStress, readTraj
export getNbStepEnergies, getNbStepStress, getNbStepAtomsFtraj
export writeEnergies, writeStress, writeFtraj
export buildingDataBase

# Read input
# Reads the input file of a CPMD simuation
#==============================================================================#
function readInputTimestep( file_input_path::T1 ) where { T1 <: AbstractString }
    # Check existence
    #-----------------------------------
    if ! isfile( file_input_path )
        print("No input file at ",file_input_path,"\n")
        return false
    end
    #-----------------------------------

    # Read input
    #-----------------------------------
    lines = utils.getLines( file_input_path )
    #-----------------------------------

    # Extract timestep
    #-----------------------------------
    timestep=0
    nb_lines=size(lines)[1]
    for line_nb=1:nb_lines
        keywords=split(lines[line_nb])
        if size(keywords)[1] > 0
            if keywords[1] == "TIMESTEP"
                timestep=parse(Float64,split(lines[line_nb+1])[1])
            end
        end
    end
    #-----------------------------------

    # Conversion to fs
    return conversion.hatime2fs*timestep
end
function readIntputStrideStress( file_input_path::T1 ) where { T1 <: AbstractString }

    # Check file existence
    #-----------------------------------
    if ! isfile( file_input_path )
        print("No input file at ",file_input_path,"\n")
        return false
    end
    #-----------------------------------

    # Read input
    #-----------------------------------
    lines = utils.getLines( file_input_path )
    #-----------------------------------

    # Extract stride of STRESS tensor
    #-----------------------------------
    stride_stress = 0
    nb_lines=size(lines)[1]
    for line_nb=1:nb_lines
        keywords=split(lines[line_nb])
        if size(keywords)[1] > 0
            if keywords[1] == "STRESS" && keywords[2] == "TENSOR"
                stride_stress = parse(Float64,split(lines[line_nb+1])[1])
            end
        end
    end
    #-----------------------------------

    return stride_stress
end
function readIntputStrideTraj( file_input_path::T1 ) where { T1 <: AbstractString }

    # Check file existence
    #-----------------------------------
    if ! isfile( file_input_path )
        print("No input file at ",file_input_path,"\n")
        return false
    end
    #-----------------------------------

    # Read input
    #-----------------------------------
    lines = utils.getLines( file_input_path )
    #-----------------------------------

    # Extract stride of STRESS tensor
    stride_traj = 0
    #-----------------------------------
    nb_lines=size(lines)[1]
    for line_nb=1:nb_lines
        keywords=split(lines[line_nb])
        if size(keywords)[1] > 0
            if keywords[1] == "TRAJECTORY"
                stride_traj = parse(Float64,split(lines[line_nb+1])[1])
            end
        end
    end
    #-----------------------------------

    return stride_traj
end
function copyInputParams( file_input::T1, handle_out::T2 ) where { T1 <: AbstractString, T2 <: IO }
    #--------------------------------------
    if ! isfile( file_input_path )
        print("No input file found at ",file_input," !\n")
        return false
    end
    #--------------------------------------

    #--------------------------------------
    file_in  = open( file_input_path )
    while ! eof( file_in )
        keywords = utils.getLineElements( file_in )
        if keywords[1] == "&ATOMS"
            break
        else
            # copying line
            #-----------------------------------------------
            utils.copyLine2file( keywords, handle_out )
            #------------------------------------------------
        end
    end
    close( file_in  )
    #--------------------------------------

    return true
end
function copyInputParams( file_input_path::T1, file_output_path::T2 ) where { T1 <: AbstractString, T2 <: AbstractString }
    file_out = open( file_out_path, "w" )
    test = copyInput( file_input_path, file_out )
    close(file_out)
    return test
end
function copyInputParams( handle_in::T1, handle_out::T2 ) where { T1 <: IO, T2 <: IO }

    while ! eof( handle_in )
        keywords = utils.getLineElements( handle_in )
        if size(keywords)[1] > 0
            if keywords[1] == "&ATOMS"
                break
            else
                # copying line
                #-----------------------------------------------
                utils.copyLine2file( keywords, handle_out )
                #------------------------------------------------
            end
        else
            write( handle_out, string("\n") )
        end
    end

    return true
end
function copyInputParams( handle_in::T1, file_output_path::T2 ) where { T1 <: IO, T2 <: AbstractString }
    file_out = open( file_output_path, "w")
    test = copyInputParams( handle_in, file_out )
    close(file_out)
    return test
end
function writeVelocities( file_out::T1, velocities::Array{T2,2} ) where { T1 <: IO, T2 <: Real }
    write( file_out, string("VELOCITIES\n") )
    nb_atoms=size(velocities)[1]
    write( file_out, string(nb_atoms," "))
    for atom=1:nb_atoms
        write( file_out, string(atom," ") )
    end
    write( file_out, string("\n") )
    for atom=1:nb_atoms
        for i=1:3
            write( file_out, string(velocities[atom,i]," ") )
        end
        write(file_out,string("\n"))
    end
    write(file_out,string("END VELOCITIES\n"))
    return true
end
function writeVelocities( file_out::T1, velocities::Array{T2,2}, nb_atoms_nb::Vector{T3} ) where { T1 <: IO, T2 <: Real, T3 <: Int }
    write( file_out, string("VELOCITIES\n") )
    write( file_out, string( size(nb_atoms_nb)[1], " " ) )
    for atom in nb_atoms_nb
        write( file_out, string( atom, " " ) )
    end
    write( file_out, string("\n") )
    for atom in nb_atoms_nb
        for i=1:3
            write( file_out, string(velocities[atom,i]," ") )
        end
        write(file_out,string("\n"))
    end
    write("END VELOCITIES\n")
    return true
end
function writePositions( file_out::T1, positions::Array{T2,2} ) where { T1 <: IO, T2 <: Real }
    nb_atoms=size(positions)[1]
    for atom=1:nb_atoms
        for i=1:3
            write(file_out,string( positions[atom,i]," " ))
        end
        write(file_out,string("\n"))
    end
    return true
end
function computeTimestep( input_path::T1 ) where { T1 <: AbstractString }
    timestep = readInputTimestep( input_path )
    if timestep == false
        return false
    end
    stride_traj = readIntputStrideTraj( input_path )
    if stride_traj == false
        return false
    end
    return stride_traj*timestep
end
#==============================================================================#

# Read Output files
#==============================================================================#
# Reading ENERGIES file
#-----------------------------------------
# Contains: Temperature, Potential Energy, Total Energy, MSD and Computing time for each step
# Structure:
# 1 line per step, per column:
# time, temperature, potential energy, total energy, MSD, Computing time
#-----------------------------------------
col_time = 1
col_temp = 3
col_poten = 4
col_entot = 5
col_msd   = 7
col_comp =  8
function getNbStepEnergies( file_path::T1 ) where { T1 <: AbstractString }
    if ! isfile( file_path )
        print("No ENERGIES file at ",file_path,"\n")
        return false
    end
    nb_step=0
    file_in = open( file_path )
    while ( ! eof(file_in) )
        temp=readline(file_in)
        nb_step += 1
    end
    return nb_step
end
function readEnergies( file_path::T1 ) where { T1 <: AbstractString }

    # Check file
    #----------------------------------------
    nb_step = getNbStepEnergies( file_path )
    if nb_step == false
        return false, false, false, false, false
    end
    #----------------------------------------

    # Array Init
    #----------------------------------------
    temp = zeros( nb_step )
    epot = zeros( nb_step )
    etot = zeros( nb_step )
    msd  = zeros( nb_step )
    comp = zeros( nb_step )
    #----------------------------------------

    # Getting data from lines
    #----------------------------------------------
    file_in = open(file_path )
    for step=1:nb_step
        line=split( readline(file_in) )
        temp[step]=parse(Float64,line[col_temp])
        epot[step]=parse(Float64,line[col_poten])
        etot[step]=parse(Float64,line[col_entot])
        msd[step]=parse(Float64,line[col_msd])
        comp[step]=parse(Float64,line[col_comp])
    end
    close(file_in)
    #----------------------------------------------

    return  temp, epot, etot, msd, comp
end
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
#------------------------------------------------------------------------------#
# Reading STRESS file
#----------------------------------------------
# Contains: Stress tensor for each step (with a possible stride)
# Structure:
# 4 lines per step
# line 1: Comment (indicates step number)
# line 2-4: stress tensor in matrix form
# Sxx Sxy Sxz
# Syx Syy Syz
# Szx Szy Szz
#----------------------------------------------
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
function readStress( file_path::T1 ) where { T1 <: AbstractString, T2 <: Int }

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
                if keywords[1] != "<<<<<<"
                    for i=1:3
                        positions[count_step,atom,i]  = parse(Float64, keywords[ i + col_start_position ] )
                        velocities[count_step,atom,i] = parse(Float64, keywords[ i + col_start_velocity ] )
                        forces[count_step,atom,i]     = parse(Float64, keywords[ i + col_start_force ] )
                    end
                end
            end
            count_step += 1
        end
    end
    close( file_in )
    #---------------------------------------------------------------------------

    return positions, velocities, forces
end
function readFttraj( file_path::T1, stride_::T2, nb_ignored::T3 ) where { T1 <: AbstractString, T2 <: Int, T3 <: Int }

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
                if keywords[1] != "<<<<<<"
                    for i=1:3
                        positions[count_step,atom,i]  = parse(Float64, keywords[ i + col_start_position ] )
                        velocities[count_step,atom,i] = parse(Float64, keywords[ i + col_start_velocity ] )
                        forces[count_step,atom,i]     = parse(Float64, keywords[ i + col_start_force ] )
                    end
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
                if keywords[1] != "<<<<<<"
                    for i=1:3
                        positions[count_step,atom,i]  = parse(Float64, keywords[ i + col_start_position ] )
                        velocities[count_step,atom,i] = parse(Float64, keywords[ i + col_start_velocity ] )
                        forces[count_step,atom,i]     = parse(Float64, keywords[ i + col_start_force ] )
                    end
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
    if writeRelaunchPositions( folder_in_target, traj[1], file_out_positions_suffix ) && writeRelaunchVelocities( folder_in_target, velocities, file_out_velocities )
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
        for atom=1:nb_atoms
            readline( file_in )
        end
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
    positions, velocities, forces = readFtraj( string( folder_in_target, "FTRAJECTORY" ) )
    if positions == false
        return false
    end
    forces = []
    nb_step = size( positions )[1]
    positions  = positions[nb_step,:,:]
    velocities = velocities[nb_step,:,:]
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
        if keywords[1] == "&END" || keywords == "VELOCITIES"
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
        writePositions( file_out, positions )
        # Ignore the atoms positions
        for atom=1:nb_atoms
            readline( file_in )
        end
    end
    write( file_out, string("\n") )
    writeVelocities( file_out, velocities )
    write( file_out, string("&END") )
    close( file_in  )
    close( file_out )
end
#-------------------------------------------------------------------------------

end
