module cpmd

using conversion
using utils

export readInputTimestep, readIntputStrideStress, readIntputStrideTraj
export readEnergy, readPressure, readStress, readTRAJ

# Read input
# Reads the input file of a CPMD simuation
#==============================================================================#
function readInputTimestep( file_input_path::T1 ) where { T1 <: AbstractString }
    # Read input
    file_in=open(file_input_path)
    lines=readlines(file_in)
    close(file_in)

    # Extract timestep
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

    # Conversion to fs
    return conversion.hatime2fs*timestep
end
function readIntputStrideStress( file_input_path::T1 ) where { T1 <: AbstractString }
    # Readinput
    file_in=open(file_input_path)
    lines=readlines(file_in)
    close(file_in)

    # Extract stride of STRESS tensor
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

    return stride_stress
end
function readIntputStrideTraj( file_input_path::T1 ) where { T1 <: AbstractString }
    # Readinput
    file_in=open(file_input_path)
    lines=readlines(file_in)
    close(file_in)

    # Extract stride of STRESS tensor
    stride_traj = 0
    nb_lines=size(lines)[1]
    for line_nb=1:nb_lines
        keywords=split(lines[line_nb])
        if size(keywords)[1] > 0
            if keywords[1] == "TRAJECTORY"
                stride_traj = parse(Float64,split(lines[line_nb+1])[1])
            end
        end
    end

    return stride_traj
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
function getEnergiesNbStep( file_path::T1 ) where { T1 <: AbstractString }
    if ! isfile( file_path )
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
function readEnergiesFile( file_path::T1 ) where { T1 <: AbstractString }
    # Check file
    if ! isfile(file_path)
        return false, false, false, false, false
    end

    # Array Init
    #----------------------------------------
    nb_steps = getEnergiesNbStep( file_path )
    temp=Vector{Real}(undef,nb_steps)
    epot=Vector{Real}(undef,nb_steps)
    etot=Vector{Real}(undef,nb_steps)
    msd=Vector{Real}(undef,nb_steps)
    comp=Vector{Real}(undef,nb_steps)
    #----------------------------------------

    # Getting data from lines
    #----------------------------------------------
    file_in = open(file_path)
    for i=1:nb_steps
        line=split( readline(file_in) )
        temp[i]=parse(Float64,line[3])
        epot[i]=parse(Float64,line[4])
        etot[i]=parse(Float64,line[5])
        msd[i]=parse(Float64,line[7])
        comp[i]=parse(Float64,line[8])
    end
    close(file_in)
    #----------------------------------------------

    return  temperature, epot, etot, msd, comp
end
function readEnergiesFile( file_path::T1, stride_::T2 ) where { T1 <: AbstractString, T2 <: Int  }
    # Check file
    if ! isfile(file_path)
        return false, false, false, false, false
    end

    # Array Init
    #----------------------------------------
    nb_steps_origin = getEnergiesNbStep( file_path )
    nb_steps = trunc( Int, nb_steps_origin/stride_ ) + 1
    temp=Vector{Real}(undef,nb_steps)
    epot=Vector{Real}(undef,nb_steps)
    etot=Vector{Real}(undef,nb_steps)
    msd=Vector{Real}(undef,nb_steps)
    comp=Vector{Real}(undef,nb_steps)
    #----------------------------------------

    # Getting data from lines
    #----------------------------------------------
    file_in = open(file_path)
    count_=1
    for step=1:nb_steps_origin
        line=split( readline(file_in) )
        if step % stride_ == 0
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
function readEnergiesFile( file_path::T1, stride_::T2, nb_ignore::T3 ) where { T1 <: AbstractString, T2 <: Int, T3 <: Int  }
    # Check file
    if ! isfile(file_path)
        return false, false, false, false, false
    end

    # Array Init
    #----------------------------------------
    nb_steps_origin = getEnergiesNbStep( file_path )
    nb_steps = trunc( Int, nb_steps_origin/stride_ ) + 1
    temp=Vector{Real}(undef,nb_steps)
    epot=Vector{Real}(undef,nb_steps)
    etot=Vector{Real}(undef,nb_steps)
    msd=Vector{Real}(undef,nb_steps)
    comp=Vector{Real}(undef,nb_steps)
    #----------------------------------------

    # Getting data from lines
    #----------------------------------------------
    file_in = open(file_path)
    count_=1
    for step=1:nb_ignore
        temp=readline(file_in)
    end
    for step=nb_ignore+1:nb_steps_origin
        line=split( readline(file_in) )
        if step % stride_ == 0 && step
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
function getNbStepStress( file_name::T1 ) where { T1 <: AbstractString }
    # Check if file exists
    if ! isfile( file_name )
        print("STRESS file does not exists!\n")
        return false
    end
    # Counting the Number of lines
    nb_line=0
    file_in = open( file_path )
    while ( ! eof(file_in) )
        temp=readline(file_in)
        nb_line += 1
    end
    # If the number of lines is not nb_line*block_size, the file is likely corrupted
    if nb_line % stress_block_size != 0
        print("STRESS file likely corrupted!\n")
        return false
    end
    # Returns number of blocks
    return Int(nb_line/stress_block_size)
end
stress_dim=3
function readStress( file_name::T1 ) where { T1 <: AbstractString, T2 <: Int }

    # Checking file exists
    if ! isfile(file_name)
        false
    end

    # Init data files
    #------------------------------------------
    nb_step=getNbStress( file_name )
    stress=zeros(Real,nb_stress_points,stress_dim,stress_dim)
    #------------------------------------------

    offset=0
    for step=1:nb_stress_points
        for i=1:3
            keywords=split(lines[1+4*(step-1)+i+offset])
            if keywords[1] == "TOTAL"
                if offset == 1
                    print("DOUBLE SIM SPOTTED at step : ",step,"\n")
                    print("LINE: ",1+4*(step-1)+i+offset,"\n")
                    return zeros(1,1), false
                end
                offset += 1
            else
                for j=1:3
                    stress[step,i,j] = parse(Float64,keywords[j])
                end
            end
        end
    end

    return stress, true
end
function readStress_( file_name::T1 ) where { T1 <: AbstractString, T2 <: Int }

    # Checking file exists
    if ! isfile(file_name)
        return zeros(1,1), false
    end
    # Reading file
    file=open(file_name);
    lines=readlines(file);
    close(file);

    # Le fichier est composés de blocs de 4 lignes:
    # 1 ligne pour le time step
    # 3 lignes pour le stress tensor
    nb_stress_points=Int(trunc(size(lines)[1]/4))
    stress=zeros(Real,nb_stress_points,3,3)
    offset=0
    for step=1:nb_stress_points
        for i=1:3
            keywords=split(lines[1+4*(step-1)+i+offset])
            if keywords[1] == "TOTAL"
                if offset == 1
                    print("DOUBLE SIM SPOTTED at step : ",step,"\n")
                    print("LINE: ",1+4*(step-1)+i+offset,"\n")
                    return zeros(1,1), false
                end
                offset += 1
            else
                for j=1:3
                    stress[step,i,j] = parse(Float64,keywords[j])
                end
            end
        end
    end

    return stress, true
end
function readStress( file_name::T1, stride::T2 ) where { T1 <: AbstractString, T2 <: Int }

    # Checking file exists
    if ! isfile(file_name)
        return zeros(1,1), false
    end

    # Reading file
    file=open(file_name);
    lines=readlines(file);
    close(file);

    # Creation de variables
    nb_stress_points=Int(trunc(size(lines)[1]/(4*stride)))
    stress=Array{Real}(nb_stress_points,3,3)
    offset=0
    # loop
    for step=1:nb_stress_points
        for i=1:3
            for j=1:3
                keywords=split(lines[1+4*(step-1)*stride+i+offset])
                if keywords[1] == "TOTAL"
                    if offset == 1
                        print("DOUBLE SIM SPOTTED at : ",step,"\n")
                        return zeros(1,1), false
                    end
                    offset+=1
                end
                stress[step,i,j] = parse(Float64,keywords[j])
            end
        end
    end
    #----------------------------------------------------------

    return stress,true
end
# Read FTRAJECTORY file
#-------------------------------------------------------------------------------
# Contains: positions, velocity and forces in atomic units for each step
# Structure:
# 1 line per atom per step
# Per line:
# atom_number x y z vx vy vz fx fy fz
#-------------------------------------
# Positions in Bohr
# Velocities in Bohr/tHart
# Forces in Ha/Bohr
#-------------------------------------
function getNbStepAtomsFTRAJ( lines::Vector{T1} ) where { T1 <: AbstractString }
    nb_lines=size(lines)[1]
    nb_restart_lines=0

    nb_atoms=0
    while split(lines[nb_atoms+1])[1] == "1"
        nb_atoms += 1
    end

    check_new=string("<<<<<<")
    nb_step=nb_lines
    for line=1:nb_lines
        if split(lines[line])[1] == check_new
            nb_step -= 1
        end
    end

    return Int(nb_step/nb_atoms), nb_atoms
end
function getNbStepAtomsFTRAJ( file::T1 ) where { T1 <: AbstractString}
    # Read file
    file_in=open(file)
    lines=readlines(file_in)
    close(file_in)

    return getNbStepAtomsFTRAJ( lines )
end
function readFTRAJ( file_input::T1 ) where { T1 <: AbstractString }

    # cols (both files):
    #   0:   natoms x nfi (natoms x 1, natoms x 2, ...)
    #   1-3: x,y,z cartesian coords [Bohr]
    #   4-6: x,y,z cartesian velocites [Bohr / thart ]
    #        thart = Hartree time =  0.024189 fs
    # FTRAJECTORY extra:
    #   7-9: x,y,z cartesian forces [Ha / Bohr]

    file_in=open( file_input )
    lines=readlines( file_in )
    close(file_in)

    nb_step, nb_atoms = getNbStepAtomsFTRAJ( lines )

    nb_lines=size(lines)[1]
    # When restarting run, a misc line:
    # <<<<<<  NEW DATA  >>>>>>
    check_new=string("<<<<<<")

    positions=zeros(nb_step,nb_atoms,3)
    velocity=zeros(nb_step,nb_atoms,3)
    forces=zeros(nb_step,nb_atoms,3)

    offset=0
    for step=1:nb_step
        for atom=1:nb_atoms
            line_element=split(lines[(step-1)*nb_atoms+atom+offset])
            if line_element[1] == check_new
                offset += 1
            else
                for i=1:3
                    positions[step,atom,i] = parse(Float64, line_element[i+1] )
                    velocity[step,atom,i]  = parse(Float64, line_element[i+4] )
                    forces[step,atom,i]    = parse(Float64, line_element[i+7] )
                end
            end
        end
    end

    return positions,velocity,forces
end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
function buildingDataBase( folder_input::T1, file_stress::T2, file_traj::T3, timestep_target::T4 ) where { T1 <: AbstractString, T2 <: AbstractString, T3 <: AbstractString, T4 <: Real }

    # Reading input file
    #---------------------------------------------------------------------------
    # Getting the stride for the STRESS file
    stride_stress = cpmd.readIntputStrideStress( file_input )
    # Getting the stride for the TRAJEC and FTRAJECTORY files
    stride_traj   = cpmd.readIntputStrideTraj(file_input )
    # Getting the timestep of the simulation
    timestep_sim  =  cpmd.readInputTimestep( file_input )
    #---------------------------------------------------------------------------

    # Computing
    #---------------------------------------------------------------------------
    n_stress = round( Int, timestep_target/( timestep_sim*stride_stress ) )
    n_traj   = round( Int, timestep_target/( timestep_sim*stride_traj ) )
    n_energy = round( Int, timestep_target/( timestep_sim*stride_traj ) )
    n_ftraj  = round( Int, timestep_target/( timestep_sim*stride_traj ) )
    #---------------------------------------------------------------------------

    # Read ENERGIES File
    #---------------------------------------------------------------------------
    file_energy=string(folder_target,"ENERGIES")
    temperature, e_pot, e_tot, msd, comp_time = readEnergyFile( file_energy )
    if ! test
        return false
    end
    #---------------------------------------------------------------------------

    # Treating ENERGIES data
    #---------------------------------------------------------------------------
    size_data_base=size(temperature)[1]
    temperature=temperature[1:n_energy:size_data_base]
    e_pot=e_pot[1:n_energy:size_data_base]
    e_tot=e_tot[1:n_energy:size_data_base]
    msd=msd[1:n_energy:size_data_base]
    comp_time=comp_time[1:n_energy:size_data_base]
    #---------------------------------------------------------------------------

    # READ STRESS file
    #---------------------------------------------------------------------------
    stress,test=cpmd.readStress(file_stress)
    if ! test
        return false
    end
    #---------------------------------------------------------------------------

    # Treating STRESS data -> Compute P and Stride
    #---------------------------------------------------------------------------
    pressure=press_stress.computePressure(stress)
    size_pressure=size(pressure)[1]
    pressure=pressure[1:n_stress:size_pressure]
    #---------------------------------------------------------------------------

    # Reading TRAJEC.xyz file
    #---------------------------------------------------------------------------
    file_traj=string(folder_in,"TRAJEC.xyz")
    traj,test=filexyz.readFastFile(file_traj)
    if ! test
        return false
    end
    #--------------------------------------------------------------------------

    # Treating TRAJ data
    #--------------------------------------------------------------------------
    size_traj=size(traj)[0]
    traj = traj[1:n_traj:size_traj]
    #--------------------------------------------------------------------------

    return temperature, e_potential, e_total, msd, comp_time, pressure, traj
end
function buildingDataBase( folder_sim::T1, timestep_target::T2 ) where { T1 <: AbstractString, T2 <: Real }

    # Determining target files paths
    #---------------------------------------------------------------------------
    file_input=string(folder_sim,"input") # Input - to get timestep + strides
    file_stress=string(folder_sim,"STRESS")  # STRESS: contains the stress tensor
    file_traj=string(folder_sim,"TRAJEC.xyz") # TRAJEC.xyz: MD Trajectory
    #---------------------------------------------------------------------------

    return buildingDataBase( file_input, file_stress, file_traj, timestep_target )
end
#-------------------------------------------------------------------------------

end
