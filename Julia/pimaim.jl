module pimaim

using LinearAlgebra

using conversion
using atom_mod
using cell_mod
using filexyz
using pdb
using periodicTable
using utils

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
    matrix2=copy(matrix)
    for i=1:3
        matrix2[i,:] /= LinearAlgebra.norm(matrix2[i,:])
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
function readCellOriginParams( path_file_len::T1, path_file_angles::T2 ) where { T1 <: AbstractString, T2 <: AbstractString }

    #-------------------------------------------------
    nb_lines_angles = utils.getNbLines( path_file_angles )
    if nb_lines_angles == false
        return false, false
    end
    nb_lines_lengths = utils.getNbLines( path_file_len )
    if nb_lines_lengths == false
        return false, false
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
    for line=1:nb_step
        key_len = split( readline( handle_in_len ) )
        key_ang = split( readline( handle_in_ang ) )
        for i=1:3
            lengths[line,i] = parse( Float64, key_len[1+i] )*conversion.bohr2Ang
            angles[line,i]  = parse( Float64, key_ang[1+i] )*tau
        end
        angles[line,:]=circshift(angles[line,:],-1)
    end
    close( handle_in_len )
    close( handle_in_ang )
    #-------------------------------------------------

    return lengths, angles
end
function traj2pdb( input_path::T1, position_path::T2, cell_angles_path::T3, cell_length_path::T4, out_path::T5 ) where { T1 <: AbstractString, T2 <: AbstractString, T3 <: AbstractString, T4 <: AbstractString, T5 <: AbstractString }
    species, species_nb = pimaim.getSpeciesAndNumber( input_path )
    traj = readPosCarTraj( position_path, species, species_nb )
    lengths, angles = readCellParams( cell_length_path, cell_angles_path )
    cells = cell_mod.makeCells( lengths, angles )
    pdb.writePdb( out_path , traj, cells )
    return true
end
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

function computeABCfromXYZ( lengths::Vector{T1}, angles::Vector{T2} ) where { T1 <: Real, T2 <: Real }
    K1 = ( cos( angles[1] ) - cos( angles[2] )*cos( angles[3] ) )/sin( angles[3] )
    K1_2 = K1*K1
    K2 = sqrt( 1 + 2*cos( angles[1] )*cos( angles[2] )*cos( angles[3] ) - cos( angles[1] )^2 - cos( angles[2] )^2 - cos( angles[3] )^2 )/sin( angles[3] )
    c = sqrt( lengths[3]/K2 )
    b = sqrt( ( lengths[2] - c*c*K1_2 )/( sin(angles[3])^2 ) )
    a = sqrt( lengths[1]^2 - b*b*cos(angles[3])*cos(angles[3]) - c*c*sin(angles[2]*angles[2]) )
    return [a,b,c]
end
function computeABCfromXYZ( lengths::Array{T1,2}, angles::Array{T2,2} ) where { T1 <: Real, T2 <: Real }
    nb_step=size(lengths)[1]
    params = zeros( nb_step )
    for step=1:nb_step
        params[step,:] = computeABCfromXYZ( lengths[step,:], angles[step,:] )
    end
    return params
end
function readCells( path_len::T1, path_angles::T2 ) where { T1 <: AbstractString, T2 <: AbstractString }
    lengths, angles = readCellParams( path_len, path_angles )
    nb_step=size(lengths)[1]
    cells = Vector{ Cell_param }(undef,nb_step)
    for step=1:nb_step
        params_lengths=computeABCfromXYZ( lengths[step,:], angles[step,:] )
        cells[step] = Cell_param( params_lengths, angles )
    end

    return cells
end

end
