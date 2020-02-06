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
            write( handle_out, string( positions_[atom,i]*conversion.ang2Bohr," " ) )
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
            if abs(matrix[j,i]) < 10^(-3)
                matrix[j,i] = 0
            end
            write( handle_out, string( matrix[j,i], " ") )
        end
        write( handle_out, string("\n") )
    end
    for i=1:3
        write(handle_out,string(lengths[i]*conversion.ang2Bohr,"\n"))
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
            if abs(matrix[j,i]) <  10^(-3)
                write( handle_in, string( 0 , " " ) )
            else
                write( handle_in, string( matrix[j,i], " " ) )
            end
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
    for atom=1:nb_atoms
        keys = split( readline( handle_in ) )
        for i=1:3
            atoms.positions[atom,i] = parse( Float64, keys[i] )
        end
        atoms.names[atom] = names_[atom]
        atoms.index[atom] = atom
    end
    #---------------------------------------------------------------------
    close( handle_in )

    #---------------------------------------------------------------------
    matrix2 = zeros(Real,3,3)
    for i=1:3
        matrix2[i,i] = LinearAlgebra.norm(matrix[i,:])
    end
    cell2=cell_mod.Cell_matrix(matrix2)
    #---------------------------------------------------------------------

    return atoms, cell2
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
function readCellParams( path_file_len::T1, path_file_angles::T2 ) where { T1 <: AbstractString, T2 <: AbstractString }

    #-------------------------------------------------
    nb_lines = utils.getNbLines( path_file_angles )
    if nb_lines == false
        return false, false
    end
    nb_lines2 = utils.getNbLines( path_file_len )
    if nb_lines2 == false
        return false, false
    end
    #-------------------------------------------------

    #-------------------------------------------------
    lengths=zeros(Real,nb_lines,3)
    angles=zeros(Real,nb_lines,3)
    #-------------------------------------------------

    #-------------------------------------------------
    tau=180/pi
    handle_in_len = open( path_file_len )
    handle_in_ang = open( path_file_angles )
    for line=1:nb_lines
        key_len = split( readline( handle_in_len ) )
        key_ang = split( readline( handle_in_ang ) )
        for i=1:3
            lengths[line,i] = parse( Float64, key_len[1+i] )*conversion.bohr2Ang
            angles[line,i]  = parse( Float64, key_ang[1+i] )*tau
        end
    end
    close( handle_in_len )
    close( handle_in_ang )
    #-------------------------------------------------

    return lengths, angles
end
function traj2pdb( input_path::T1, position_path::T2, cell_angles_path::T3, cell_length_path::T4, out_path::T5 ) where { T1 <: AbstractString, T2 <: AbstractString, T3 <: AbstractString, T4 <: AbstractString, T5 <: AbstractString }
    species, species_nb = pimaim.getSpeciesAndNumber( input_path )
    nb_atoms = sum(species_nb)
    names_atoms = atom_mod.buildNames( species, species_nb )
    positions=readPositions( position_path, nb_atoms )
    index_atoms=[1:nb_atoms;]
    nb_step = size(positions)[1]
    lengths, angles = readCellParams( cell_length_path, cell_angles_path )
    traj = atom_mod.makeTrajAtomList( positions,  names_atoms, index_atoms )
    cells = cell_mod.makeCells( lengths, angles )
    pdb.writePdb( out_path , traj, cells )
    return true
end

end
