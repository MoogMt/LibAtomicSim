module pimaim

using conversion
using LinearAlgebra
using atom_mod
using cell_mod
using filexyz
using pdb
using periodicTable
using utils

# CONTAINS
# - Readers:
# getSpeciesAndNumber( path_file::AbstractString ) : Reads runtime.inpt and extracts species and their number of atoms
# readPositions( path_file::AbstractString, nb_atoms::Int ) : Old - to be deleted
# readPositionsUpToCrash( path_file::AbstractString, nb_atoms ) : Old to be deleted
# readPositions( path_file::AbstractString, species::Vector(String), nb_species::Vector{Int}, cells::Vector{Cell_params} )

# Reads runtime.inpt to get the species and number of species
#-------------------------------------------------------------------------------
function getSpeciesAndNumber( path_file::T1 ) where { T1 <: AbstractString }
    if ! isfile( path_file )
        return false, false
    end
    handle_in = open( path_file )
    utils.skipLines( handle_in, 4 )
    nb_species=parse( Int, split( readline( handle_in ) )[1] )
    #---------------------------------------
    species=Vector{AbstractString}(undef,nb_species)
    species_line = split( readline( handle_in ), "," )
    for i_spec = 1:nb_species-1
        species[i_spec] = species_line[i_spec]
    end
    species[ nb_species ] = split(species_line[nb_species])[1]  # Avoid pesky comment
    #---------------------------------------
    species_nb = zeros(Int, nb_species )
    species_nb_line = split( readline( handle_in ),"," )
    for i_spec = 1:nb_species-1
        species_nb[i_spec] = parse(Int, species_nb_line[i_spec] )
    end
    species_nb[ nb_species ] = parse(Int, split(species_nb_line[nb_species])[1] )  # Avoid pesky comment
    #---------------------------------------
    close( handle_in )
    return species, species_nb
end
#-------------------------------------------------------------------------------
function readPositions( path_file::T1, nb_atoms::T2 ) where { T1 <: AbstractString, T2 <: Int }

    nb_lines = utils.getNbLines( path_file )
    if nb_lines == false
        return false
    end

    if nb_lines % nb_atoms != 0
        print("The number of lines is not a multiple of the number of atoms, either the number of atoms is wrong or the file is corrupted.\n")
        print("File: ",path_file,"\n")
        return false
    end

    nb_step = Int(nb_lines/nb_atoms)
    positions=zeros(Real,nb_step,nb_atoms,3)

    handle_in = open( path_file )
    for step=1:nb_step
        for atom=1:nb_atoms
            keys = split( readline( handle_in ) )
            for i=1:3
                positions[step,atom,i] = parse( Float64, keys[i] )*conversion.bohr2Ang
            end
        end
    end
    close( handle_in )

    return positions
end
function readPositionsUpToCrash( path_file::T1, nb_atoms::T2 ) where { T1 <: AbstractString, T2 <: Int }

    nb_lines = utils.getNbLines( path_file )
    if nb_lines == false
        return false
    end

    nb_step = Int( trunc( nb_lines/nb_atoms ) )
    positions=zeros(Real,nb_step,nb_atoms,3)

    handle_in = open( path_file )
    for step=1:nb_step
        for atom=1:nb_atoms
            keys = split( readline( handle_in ) )
            for i=1:3
                positions[step,atom,i] = parse( Float64, keys[i] )*conversion.bohr2Ang
            end
        end
    end
    close( handle_in )

    return positions
end
function readPosCar( path_file::T1 ) where { T1 <: AbstractString }

    #---------------------------------------------------------------------
    nb_lines = utils.getNbLines( path_file )
    if nb_lines == false
        return false, false
    end
    #---------------------------------------------------------------------

    handle_in = open( path_file )
    utils.skipLines( handle_in, 2) # Skip the three first lines
    #---------------------------------------------------------------------
    matrix=zeros(Real,3,3)
    for i=1:3
        keyword=split( readline( handle_in) )
        for j=1:3
            matrix[i,j] = parse(Float64,keyword[j])
        end
    end
    cell = cell_mod.Cell_matrix(matrix)
    #---------------------------------------------------------------------
    matrix2  = copy( matrix )
    for i=1:3
        matrix2[i,:] /= LinearAlgebra.norm(matrix2[i,:])
    end
    if matrix2[2,1] == 0
        matrix2 = transpose( matrix2 )
    end
    #---------------------------------------------------------------------
    species=split( readline( handle_in ) )
    nb_species_ = split( readline(handle_in ) )
    nb_species = zeros( Int, size(species)[1] )
    for i_spec=1:size(species)[1]
        nb_species[i_spec] = parse( Int, nb_species_[i_spec] )
    end
    names_ = atom_mod.buildNames( species, nb_species )
    #---------------------------------------------------------------------
    readline( handle_in ) # skip
    nb_atoms=sum(nb_species)
    atoms=AtomList(nb_atoms)
    temp = zeros(Real,3)
    for atom=1:nb_atoms
        keys = split( readline( handle_in ) )
        for i=1:3
            temp[i] = parse( Float64, keys[i] )
        end
        for i=1:3
            for j=1:3
                atoms.positions[atom,i] += matrix2[j,i]*temp[j]
            end
        end
        atoms.names[atom] = names_[atom]
        atoms.index[atom] = atom
    end
    #---------------------------------------------------------------------
    close( handle_in )

    cell.matrix=transpose(cell.matrix)

    return atoms, cell
end
function readPosCarTraj( path_file::T1, species::Vector{T2}, nb_species::Vector{T3} ) where { T1 <: AbstractString, T2 <: AbstractString, T3 <: Int }

    #---------------------------------------------------------------------
    nb_lines = utils.getNbLines( path_file )
    if nb_lines == false
        return false
    end
    #---------------------------------------------------------------------

    #---------------------------------------------------------------------
    nb_atoms = sum(nb_species)
    names_=atom_mod.buildNames( species, nb_species )
    if nb_lines % nb_atoms != 0
        print("Problem within file: ",path_file," :\n")
        print("Number of lines is not a multiple of the number of atoms given.\n")
        return false
    end
    nb_step = Int(nb_lines/nb_atoms)
    traj = Vector{AtomList}( undef, nb_step )
    #---------------------------------------------------------------------

    #---------------------------------------------------------------------
    handle_in = open( path_file )
    for step=1:nb_step
        traj[step] = atom_mod.AtomList(nb_atoms)
        for atom=1:nb_atoms
            keys = split( readline( handle_in ) )
            for i=1:3
                traj[step].positions[atom,i] = parse( Float64, keys[i] )*conversion.bohr2Ang
            end
            traj[step].names[atom] = names_[atom]
            traj[step].index[atom] = atom
        end
    end
    #---------------------------------------------------------------------
    close( handle_in )

    return traj
end
function readTrajPosCar( input_path::T1, poscar_path::T2, cell_length_path::T3, cell_angles_path::T4 ) where { T1 <: AbstractString, T2 <: AbstractString, T3 <: AbstractString, T4 <: AbstractString }
    species, species_nb = pimaim.getSpeciesAndNumber( input_path )
    if species == false || species_nb == false
        return false, false
    end
    traj = readPosCarTraj( poscar_path, species, species_nb )
    if traj == false
        return false, false
    end
    cells = pimaim.readCellParams( cell_length_path, cell_angles_path )
    if cells == false
        return false, false
    end
    return traj, cells
end
function readPositions( path_file::T1, species::Vector{T2}, nb_species::Vector{T3}, cells::Vector{T4} ) where { T1 <: AbstractString, T2 <: AbstractString, T3 <: Int, T4 <: cell_mod.Cell_param }

    #---------------------------------------------------------------------
    nb_lines = utils.getNbLines( path_file )
    if nb_lines == false
        # utils.getNbLines will have return the error message
        return false
    end
    #---------------------------------------------------------------------

    #---------------------------------------------------------------------
    nb_atoms = sum(nb_species)
    names_=atom_mod.buildNames( species, nb_species )
    if nb_lines % nb_atoms != 0
        print("Problem within file: ",path_file," :\n")
        print("Number of lines is not a multiple of the number of atoms given.\n")
        return false
    end
    nb_step = Int(nb_lines/nb_atoms)
    traj = Vector{AtomList}( undef, nb_step )
    #---------------------------------------------------------------------

    #---------------------------------------------------------------------
    handle_in = open( path_file )
    for step=1:nb_step
        box_len=zeros(Real,3)
        traj[step] = atom_mod.AtomList(nb_atoms)
        matrix2 = copy( cell_mod.params2Matrix( cells[step] ).matrix )
        for i=1:3
            box_len[i] = LinearAlgebra.norm( matrix2[:,i] )*conversion.ang2Bohr
        end
        for atom=1:nb_atoms
            temp = zeros(Real,3)
            keys = split( readline( handle_in ) )
            for i=1:3
                temp[i] = parse(Float64,keys[i])/box_len[i]
            end
            for i=1:3
                for j=1:3
                    traj[step].positions[atom,i] += temp[j]*matrix2[i,j]
                end
            end
            traj[step].names[atom] = names_[atom]
            traj[step].index[atom] = atom
        end
    end
    #---------------------------------------------------------------------
    close( handle_in )

    return traj
end
function readPositionsUpToCrash( path_file::T1, species::Vector{T2}, nb_species::Vector{T3}, cells::Vector{T4} ) where { T1 <: AbstractString, T2 <: AbstractString, T3 <: Int, T4 <: cell_mod.Cell_param }

    #---------------------------------------------------------------------
    nb_lines = utils.getNbLines( path_file )
    if nb_lines == false
        return false
    end
    #---------------------------------------------------------------------

    #---------------------------------------------------------------------
    nb_atoms = sum(nb_species)
    names_=atom_mod.buildNames( species, nb_species )
    nb_step = Int( trunc( nb_lines/nb_atoms ) )-1
    if nb_step < 1
        return false
    end
    if size(cells)[1] < nb_step
        nb_step = size(cells)[1]
    end
    traj = Vector{AtomList}( undef, nb_step )
    #---------------------------------------------------------------------

    #---------------------------------------------------------------------
    handle_in = open( path_file )
    for step=1:nb_step
        box_len=zeros(Real,3)
        traj[step] = atom_mod.AtomList(nb_atoms)
        matrix2 = copy( cell_mod.params2Matrix( cells[step] ).matrix )
        for i=1:3
            box_len[i] = LinearAlgebra.norm( matrix2[:,i] )*conversion.ang2Bohr
        end
        for atom=1:nb_atoms
            temp = zeros(Real,3)
            keys = split( readline( handle_in ) )
            for i=1:3
                temp[i] = parse(Float64,keys[i])/box_len[i]
            end
            for i=1:3
                for j=1:3
                    traj[step].positions[atom,i] += temp[j]*matrix2[i,j]
                end
            end
            traj[step].names[atom] = names_[atom]
            traj[step].index[atom] = atom
        end
    end
    #---------------------------------------------------------------------
    close( handle_in )

    return traj
end
function readXV( path_file::T1 ) where { T1 <: AbstractString }
    matrix=zeros(Real,3,3)
    handle_in = open( path_file )
    for i=1:3
        line = split(readline( handle_in ))
        for j=1:3
            matrix[i,j] = parse(Float64, line[j] )*conversion.bohr2Ang
        end
    end
    cell = cell_mod.cellMatrix2Params( cell_mod.Cell_matrix(matrix) )
    nb_atoms = parse(Int64, split( readline( handle_in ) )[1] )
    atoms = atom_mod.AtomList( nb_atoms )
    for atom=1:nb_atoms
        keyword = split( readline( handle_in ) )
        atoms.names[atom] = periodicTable.z2Names( parse(Int, keyword[2] ) )
        atoms.index[atom] = atom
        for i=1:3
            atoms.positions[atom,i] = parse(Float64, keyword[2+i] )*conversion.bohr2Ang
        end
    end
    close( handle_in )
    return atoms, cell
end
function readRestart( restart_path::T1, runtime_path::T2 ) where { T1 <: AbstractString, T2 <: AbstractString }

    #-----------------------------------------------------------------------
    species, nb_element_species = getSpeciesAndNumber( runtime_path )
    n_species = size(species)[1]
    nb_atoms = sum( nb_element_species )
    atoms_names = Vector{AbstractString}(undef,nb_atoms)
    for i_spec=1:n_species
        offset_specie = sum( nb_element_species[1:i_spec-1] )
        for atom_spec=1:nb_element_species[i_spec]
            atoms_names[ offset_specie + atom_spec ] = species[ i_spec ]
        end
    end
    #-----------------------------------------------------------------------

    #-----------------------------------------------------------------------
    handle_in = open( restart_path )
    check_position  = split( readline( handle_in ) )[1]
    check_velocity  = split( readline( handle_in ) )[1]
    check_forces    = split( readline( handle_in ) )[1]
    check_polarity  = split( readline( handle_in ) )[1]
    #-----------------------------------------------------------------------

    #-----------------------------------------------------------------------
    atoms = false
    if check_position == "T"
        atoms = atom_mod.AtomList( nb_atoms )
        for atom=1:nb_atoms
            keyword = split( readline( handle_in ) )
            atoms.index[atom] = atom
            atoms.names[atom] = atoms_names[atom]
            for i=1:3
                atoms.positions[atom,i] = parse(Float64, keyword[i] )*conversion.bohr2Ang
            end
        end
    end
    #-----------------------------------------------------------------------

    #-----------------------------------------------------------------------
    velocity = false
    if check_velocity == "T"
        velocities = zeros(Real, nb_atoms, 3)
        for atom=1:nb_atoms
            keyword = split( readline( handle_in ) )
            for i=1:3
                velocities[atom,i] = parse(Float64, keyword[i] )
            end
        end
    end
    #-----------------------------------------------------------------------

    #-----------------------------------------------------------------------
    forces = false
    if check_forces == "T"
        forces = zeros(Real, nb_atoms, 3)
        for atom=1:nb_atoms
            keyword = split( readline( handle_in ) )
            for i=1:3
                forces[ atom, i ] = parse(Float64, keyword[i] )
            end
        end
    end
    #-----------------------------------------------------------------------

    #-----------------------------------------------------------------------
    nb_runs = parse(Int64, split( readline( handle_in ) )[1])
    runs_nb_step = zeros( nb_runs )
    for i=1:nb_runs
        runs_nb_step[i] = parse(Int64, split( readline(handle_in) )[1] )
    end
    #-----------------------------------------------------------------------

    #-----------------------------------------------------------------------
    polarization = zeros( 30 )
    for i=1:30
        polarization[i] = parse(Float64, split( readline(handle_in) )[1] )
    end
    #-----------------------------------------------------------------------

    #-----------------------------------------------------------------------
    quad = zeros( 2, 2 )
    for i=1:2
        keys = split( readline( handle_in ) )
        for j=1:2
            quad[i,j] = parse(Float64, keys[j] )
        end
    end
    #-----------------------------------------------------------------------

    #-----------------------------------------------------------------------
    for i=1:3
        readline( handle_in )
    end
    #-----------------------------------------------------------------------

    #-----------------------------------------------------------------------
    cell_matrix = zeros(Real,3,3)
    for i=1:3
        keyword = split( readline(handle_in) )
        for j=1:3
            cell_matrix[i,j] = parse( Float64, keyword[j] )
        end
    end
    boxlens = zeros(3)
    for i=1:3
        key = split( readline( handle_in ) )
        boxlens[i] = parse(Float64, key[1] )
    end
    for i=1:3
        for j=1:3
            cell_matrix[i,j] *= cell_matrix[i,j]*boxlens[i]*conversion.bohr2Ang
        end
    end
    cell = cell_mod.Cell_matrix(cell_matrix)
    #-----------------------------------------------------------------------

    return atoms, velocity, forces, polarization, quad, cell, runs_nb_step
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
function readCellParams( path_file_len::T1, path_file_angles::T2 ) where { T1 <: AbstractString, T2 <: AbstractString }

    #-------------------------------------------------
    nb_lines_angles = utils.getNbLines( path_file_angles )
    if nb_lines_angles == false
        return false
    end
    nb_lines_lengths = utils.getNbLines( path_file_len )
    if nb_lines_lengths == false
        return false
    end
    nb_step=min( nb_lines_angles, nb_lines_lengths )
    if nb_lines_angles != nb_lines_lengths
        print("Missmatch between number of lengths and angles, using the minimum to proceed.\n")
    end
    #-------------------------------------------------

    #-------------------------------------------------
    lengths = zeros(Real, nb_step , 3 )
    angles  = zeros(Real, nb_step , 3 )
    #-------------------------------------------------

    #-------------------------------------------------
    tau=180/pi
    handle_in_len = open( path_file_len )
    handle_in_ang = open( path_file_angles )
    cells = Vector{ cell_mod.Cell_param }( undef, nb_step )
    for line=1:nb_step
        key_len = split( readline( handle_in_len ) )
        key_ang = split( readline( handle_in_ang ) )
        for i=1:3
            lengths[line,i] = parse( Float64, key_len[1+i] )*conversion.bohr2Ang
            angles[line,i]  = parse( Float64, key_ang[1+i] )*tau
        end
        stock=angles[line,1]
        angles[line,1]=angles[line,3]
        angles[line,3]=stock
        cells[ line ] = cell_mod.Cell_param( lengths[line,:], angles[line,:] )
        #angles[line,:]=circshift(angles[line,:],1)
    end
    close( handle_in_len )
    close( handle_in_ang )
    #-------------------------------------------------

    return cells
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
function readEnergy( path_file::T1 ) where { T1 <: AbstractString }

    #---------------------------------------------------------------------
    nb_lines = utils.getNbLines( path_file )
    if nb_lines == false
        return false
    end
    #---------------------------------------------------------------------

    enpot = zeros(Real,nb_lines)
    enkin = zeros(Real,nb_lines)
    handle_in = open( path_file )
    for step=1:nb_lines
        keyword = split(readline(handle_in))
        enpot[step] = parse(Float64, keyword[2] )
        enkin[step] = parse(Float64, keyword[3] )
    end
    close(handle_in)

    return enpot.+enkin, enpot, enkin
end
function readEnthalpy( path_file::T1 ) where { T1 <: AbstractString }

    #---------------------------------------------------------------------
    nb_lines = utils.getNbLines( path_file )
    if nb_lines == false
        return false
    end
    #---------------------------------------------------------------------

    enthalpy = zeros(Real,nb_lines)
    pv = zeros(Real,nb_lines)
    handle_in = open( path_file )
    for step=1:nb_lines
        keyword = split(readline(handle_in))
        pv[step] = parse(Float64, keyword[2] )
        enthalpy[step] = parse(Float64, keyword[3] )
    end
    close(handle_in)

    return enthalpy, pv
end
function readPressure( path_file::T1 ) where { T1 <: AbstractString }

    #---------------------------------------------------------------------
    nb_lines = utils.getNbLines( path_file )
    if nb_lines == false
        return false
    end
    #---------------------------------------------------------------------

    pressure = zeros(Real,nb_lines)
    handle_in = open( path_file )
    for step=1:nb_lines
        pressure[step] = parse(Float64, split(readline(handle_in))[2] )*conversion.au2Gpa
    end
    close(handle_in)

    return pressure
end
function readFullOutput( path_file::T1 ) where { T1 <: AbstractString }
    nb_lines = utils.getNbLines( path_file)
    if nb_lines == false
        return false
    end
    if nb_lines % 7 != 0
        print("Issue with fulloutput.dat at ", path_file,"!\n")
        return false
    end
    nb_step = Int( nb_lines/7 )
    handle_in = open( path_file )
    data = zeros(nb_step,6)
    for step=1:nb_step
        for i=1:6
            test=readline( handle_in )
        end
        keys = readline( handle_in )
        for i=1:6
            data[ step, i ] = parse( Float64, split(keys)[i+1] )
        end
    end
    close(handle_in)
    return data
end
function readf3( path_file::T1 ) where { T1 <: AbstractString }
    handle_in = open( path_file )
    nb_step = 0
    while ! eof( handle_in )
        readline( handle_in )
        nb_step += 1
    end
    seekstart( handle_in )
    col_nb = 9
    readline( handle_in )
    nb_step = nb_step-1
    ring_data = zeros( nb_step, col_nb )
    for step = 1:nb_step
        keys = split( readline( handle_in ) )
        for col=1:col_nb
            ring_data[ step, col ] = parse(Float64, keys[col] )
        end
    end
    close( handle_in )
    return ring_data
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
function extractTemperature( data::Array{T1,2} ) where { T1 <: Real }
    return data[:,6]
end
function extractPressure( data::Array{T1,2} ) where { T1 <: Real }
    return data[:,5]
end
function extractVolume( data::Array{T1,2} ) where { T1 <: Real }
    return data[:,4]
end
function extractTotalEnergy( data::Array{T1,2} ) where { T1 <: Real }
    return data[:,3]
end
function extractKineticEnergy( data::Array{T1,2} ) where { T1 <: Real }
    return data[:,2]
end
function extractPotentialEnergy( data::Array{T1,2} ) where { T1 <: Real }
    return data[:,1]
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
function writeFullData( handle_out::T1, data::Array{T2,2} ) where { T1 <: IO, T2 <: Real }
    nb_step = size( data )[1]
    n_dim = size( data )[2]
    for step=1:nb_step
        for i=1:n_dim
            Base.write( handle_out, string( data[step,i], " " ) )
        end
        Base.write( handle_out, string( "\n" ) )
    end
    return true
end
function writeFullData( file_out::T1, data::Array{T1,2} ) where { T1 <: AbstractString, T2 <: Real }
    handle_out = open( file_out, "w"  )
    test=writeFullData( handle_out, data )
    close( handle_out )
    return test
end
function writeRestart( path_file::T1, atoms::T2, cell::T3 ) where { T1 <: AbstractString, T2 <: atom_mod.AtomList, T3 <: cell_mod.Cell_matrix }
    nb_atoms=atom_mod.getNbAtoms(atoms)
    handle_out=open( path_file , "w" )
    write( handle_out, string("T\n") )
    write( handle_out, string("F\n") )
    write( handle_out, string("F\n") )
    write( handle_out, string("F\n") )
    positions_ = copy(atoms.positions)
    positions_ = cell_mod.getTransformedPosition( positions_, inv(cell.matrix) )
    for i=1:3
        positions_[:,i] *= LinearAlgebra.norm( cell.matrix[:,i] )
    end
    for atom=1:nb_atoms
        for i=1:3
            write( handle_out, string( round(positions_[atom,i]*conversion.ang2Bohr,digits=3)," " ) )
        end
        write( handle_out, string("\n") )
    end
    matrix=copy(cell.matrix)
    lengths=cell_mod.cellMatrix2Params(cell).length
    for i=1:3
        matrix[:,i] = matrix[:,i]/lengths[i]
    end
    for i=1:3
        for j=1:3
            write( handle_out, string( round(matrix[i,j],digits=3), " ") )
        end
        write( handle_out, string("\n") )
    end
    for i=1:3
        write(handle_out,string(round(LinearAlgebra.norm(cell.matrix[:,i])*conversion.ang2Bohr,digits=3),"\n"))
    end
    close(handle_out)
    return true
end
function writeCrystalCell( path_file::T1, atoms::T2, cell::T3 ) where { T1 <: AbstractString, T2 <: atom_mod.AtomList, T3 <: cell_mod.Cell_matrix }

    lengths=cell_mod.cellMatrix2Params(cell).length
    matrix = copy( cell.matrix )
    species = atom_mod.getSpecies( atoms )
    for i=1:3
        matrix[:,i] = matrix[:,i]/lengths[i]
    end

    handle_in = open( path_file, "w" )
    for i=1:3
        for j=1:3
            write( handle_in, string( round(matrix[j,i],digits=3), " " ) )
        end
        write( handle_in, string("\n") )
    end
    for i=1:3
        write( handle_in, string( lengths[i], "\n")  )
    end
    for i=1:3
        write( handle_in, "1\n" )
    end
    for i=1:size(species)[1]
        write( handle_in, string( periodicTable.names2Z( species[i] ), "\n" ) )
        write( handle_in, string( species[i], "_quartz.mat\n" ) )
    end
    for i=1:3
        write( handle_in, string( lengths[i], "\n")  )
    end
    close(handle_in)
    return true
end
function writeAcellTxt( path_file::T1, cells::Vector{T2} ) where { T1 <: AbstractString, T2 <: cell_mod.Cell_param }
    handle_out=open( path_file, "w" )
    nb_step = size(cells)[1]
    for step = 1:nb_step
        cell_matrix = cell_mod.params2Matrix(cells[step]).matrix
        for i=1:3
            for j=1:3
                write( handle_out, string( round(cell_matrix[i,j],digits=3), " " ) )
            end
            write( handle_out, "\n" )
        end
    end
    close(handle_out)
    return true
end
#------------------------------------------------------------------------------

end
