module pdb

export getNbSteps, readStructure, readTraj
export writePdb, writePdbPlumed

using utils
using atom_mod
using cell_mod
using periodicTable

# Reading file
#==============================================================================#
function getNbSteps( file_path::T1 ) where { T1 <: AbstractString }

    #-----------------------------------------
    if ! isfile( file_path )
        print("No pdb file found at ",file_path," !\n")
        return false
    end
    #-----------------------------------------

    #-----------------------------------------
    nb_step  = 0
    nb_step2 = 0
    file_in = open( file_path )
    while !eof(f)
        keyword1 = split(readline(f))[1]
        if keyword1 == "END"
            nb_step += 1
        elseif keyword1 == "CRYST1"
            nb_step2 +=1
        end
    end
    close( file_in )
    #-----------------------------------------

    if count_step == count_step2
        return nb_step
    else
        print("PDB file at ",file_path," is probably corrupted.\n")
        return false
    end
end
function getNbAtoms( file_path::T1 ) where { T1 <: AbstractString }

    #-----------------------------------------
    if ! isfile( file_path )
        print("No pdb file found at ",file_path," !\n")
        return false
    end
    #-----------------------------------------

    #-----------------------------------------
    nb_atoms  = 0
    file_in = open( file_path )
    while !eof(f)
        keyword1 = split(readline(f))[1]
        if keyword1 == "ATOM"
            nb_atoms += 1
        elseif keyword1 == "END"
            break
        end
    end
    close( file_in )
    #-----------------------------------------

    return nb_atoms
end
function readStructure( file_path::T1 ) where { T1 <: AbstractString }
  #--------------
  # Reading file
  #----------------------
  lines = utils.getLines( file_path )
  #------------------------

  #----------------------------------------------------
  # Reading informations about cell and number of atoms
  #----------------------------------------------------
  cell = cell_mod.Cell_param()
  nb_atoms=0
  for line in lines
    if split(line)[1] == "CRYST1"
      for i=1:3
        cell.length[i]=parse(Float64,split(line)[i+1])
        cell.angles[i]=parse(Float64,split(line)[i+4])
      end
    elseif split(line)[1] == "ATOM"
      nb_atoms+=1
    end
  end

  #----------------------------------------------------

  #---------------------------------
  # Reading atomic informations
  #---------------------------------------------------------------------
  atoms=atom_mod.AtomMolList(nb_atoms)
  count=1
  for line in lines
    if split(line)[1] == "ATOM"
      atoms.atom_names[count]=split(line)[3]
      atoms.atom_index[count]=parse(Float64,split(line)[2])
      atoms.mol_names[count]=split(line)[4]
      atoms.mol_index[count]=parse(Float64,split(line)[6])
      atoms.positions[count,:]=[ parse(Float64,split(line)[7]), parse(Float64,split(line)[8]), parse(Float64,split(line)[9]) ]
      count+=1;
    end
  end
  #---------------------------------------------------------------------

  return atoms, cell
end
function readTraj( file_path::T1 ) where { T1 <: AbstractString }
  #--------------
  # Reading file
  #----------------------
  lines = utils.getLines( file_path )
  #------------------------

  #----------------------------------------------------
  # Reading informations about cell and number of atoms
  #----------------------------------------------------
  cell = cell_mod.Cell_param()
  nb_atoms=0
  for line in lines
    if split(line)[1] == "CRYST1"
      for i=1:3
        cell.length[i]=parse(Float64,split(line)[i+1])
        cell.angles[i]=parse(Float64,split(line)[i+4])
      end
    elseif split(line)[1] == "ATOM"
      nb_atoms+=1
    end
  end

  #----------------------------------------------------

  #---------------------------------
  # Reading atomic informations
  #---------------------------------------------------------------------
  atoms=atom_mod.AtomMolList(nb_atoms)
  count=1
  for line in lines
    if split(line)[1] == "ATOM"
      atoms.atom_names[count]=split(line)[3]
      atoms.atom_index[count]=parse(Float64,split(line)[2])
      atoms.mol_names[count]=split(line)[4]
      atoms.mol_index[count]=parse(Float64,split(line)[6])
      atoms.positions[count,:]=[ parse(Float64,split(line)[7]), parse(Float64,split(line)[8]), parse(Float64,split(line)[9]) ]
      count+=1;
    end
  end
  #---------------------------------------------------------------------

  return atoms, cell
end
#==============================================================================#

# Write file
#==============================================================================#
function addJustifyRight( max_column::T1, string_target::T2, string_toadd::T3 ) where  { T1 <: Int, T2 <: AbstractString, T3 <: AbstractString }
    string_target = utils.spaces( string_target, max_column-length(string_target)-length(string_toadd) )
    return string( string_target, string_toadd)
end
function addJustifyLeft( max_column::T1, string_target::T2, string_toadd::T3 ) where { T1 <: Int, T2 <: AbstractString, T3 <: AbstractString }
    string_target = utils.spaces( string_target, max_column-length(string_target) )
    return string( string_target, string_toadd)
end
function writeCRYST1(  handle_out::IO, lengths::Vector{Real}, angles::Vector{Real},
    round_length::Int=3,
    round_angles::Int=2,
    space_group::AbstractString="P 1",
    z_number::Int=1,
    end_a::Int=15,
    end_b::Int=24,
    end_c::Int=33,
    end_alpha::Int=40,
    end_beta::Int=47,
    end_gamma::Int=54,
    start_sg::Int=56,
    end_z::Int=70 )

    # Role: write CRYST1 to file
    # Parametric, in order to adapt to the various PDB format in existence

    # Arguments
    # handle_out: IO stream to target file
    # lengths:   lengths parameters of the cell
    # angles:    angles parameters of the cell

    # Optionnal Arguments
    # round_length: number of digits for lengths rounding
    # round_angles: number of digits for angles rounding
    # space_group: Space group in Hermann Mauguin notation, by default P 1 works...
    #              never seen a case where something else is used 28/01/2020 (M Moog)
    # z_number: Number of polumeric chains or in hetero polymers
    #           the n umber of occurences of the most populous chain
    # end_a:  Maximum colum for a (justify right for a )
    # end_b;  Maximum colum for b (justify right for b )
    # end_c;  Maximum colum for c (justify right for b )
    # end_alpha;  Maximum colum for alpha (justify right for b )
    # end_beta;  Maximum colum for beta (justify right for b )
    # end_gamma;  Maximum colum for gamma (justify right for b )
    # start_sg;  Min column for Space Group
    # end_z;  Maximum colum for gamma (justify right for b )

    # lengths
    a = string( round( lengths[1], digits=round_length ) )
    b = string( round( lengths[2], digits=round_length ) )
    c = string( round( lengths[3], digits=round_length ) )
    # angles
    alpha = string( round( angles[1], digits=round_angles ) )
    beta  = string( round( angles[2], digits=round_angles ) )
    gamma = string( round( angles[3], digits=round_angles ) )
    # z number
    z_number=string(z_number)

    # Writting CRYST1 info -> Cell information
    cryst1 = string("CRYST1")
    # Cell lengths (right justified)
    cryst1 = addJustifyRight( end_a, cryst1, a )
    if cryst1 == false
        print("Error writting CRYST1 line (a).\n")
        return false
    end
    cryst1 = addJustifyRight( end_b, cryst1, b )
    if cryst1 == false
        print("Error writting CRYST1 line (b).\n")
        return false
    end
    cryst1 = addJustifyRight( end_c, cryst1, c )
    if cryst1 == false
        print("Error writting CRYST1 line (c).\n")
        return false
    end
    # Angles (right justified)
    cryst1 = addJustifyRight( end_alpha, cryst1, alpha )
    if cryst1 == false
        print("Error writting CRYST1 line (alpha).\n")
        return false
    end
    cryst1 = addJustifyRight( end_beta,  cryst1, beta  )
    if cryst1 == false
        print("Error writting CRYST1 line (beta).\n")
        return false
    end
    cryst1 = addJustifyRight( end_gamma, cryst1, gamma )
    if cryst1 == false
        print("Error writting CRYST1 line (gamma).\n")
        return false
    end
    # Space Group (left justified)
    cryst1 = addJustifyLeft( start_sg, cryst1, space_group )
    if cryst1 == false
        print("Error writting CRYST1 line (Space Group).\n")
        return false
    end
    # Z number (right justified)
    cryst1 = addJustifyRight( end_z, cryst1, z_number )
    if cryst1 == false
        print("Error writting CRYST1 line (z).\n")
        return false
    end
    cryst1=string(cryst1,"\n")

    # Writting line
    Base.write(handle_out,cryst1)

    return true
end
#------------------------------------------------------------------------------#
function writePdb( file::T1, atoms::T2, cell::T3 ) where { T1 <: AbstractString , T2 <: atom_mod.AtomList, T3 <: cell_mod.Cell_param }

    # Writes a PDB file according to standard format (2011)

    handle_out=open(file,"w")

    if ! writeCRYST1( handle_out, cell.length, cell.angles )
        print("Error writting cell information!\n")
        return false
    end


    # Atomic Positions
    nb_atoms = size(atoms.names)[1]
    for i=1:nb_atoms
        atom="ATOM"
        # Right justified
        atom=utils.spaces(atom,11-length(string(atoms.index[i]))-length(atom))
        atom=string(atom,atoms.index[i])
        # If atom name is 1 length, start on 13, otherwise 14
        # atom names are left justified here
        if length(atoms.names[i])==1
            atom=utils.spaces(atom,13-length(atom))
        else
            atom=utils.spaces(atom,14-length(atom))
        end
        atom=string(atom,atoms.names[i])

        # Alternate location indicator ( default blank )
        alt_loc=string(" ")
        atom=utils.spaces(atom,17-length(atom)-length(alt_loc))
        atom=string(atom," ")

        # Residue name (3 col max) righ justified?
        residue_name=string("   ")
        atom=utils.spaces(atom,20-length(atom)-length(residue_name))
        atom=string(atom,residue_name)

        # Chain Identifier
        # Default is "X"
        chain_id=string("X")
        atom=utils.spaces(atom,22-length(atom)-length(chain_id))
        atom=string(atom,"X")

        # Residue Sequence Nb - Molecule nb
        mol_nb=string("1")
        atom=utils.spaces(atom,26-length(atom)-length(mol_nb))
        atom=string(atom,mol_nb)

        # Code for insertion of residues
        code_insertion=string(" ")
        atom=utils.spaces(atom,27-length(atom)-length(code_insertion))
        atom=string(atom,code_insertion)

        # Positions are right justified
        # X
        x=string( round(atoms.positions[i,1], digits=3 ) )
        atom=utils.spaces( atom, 38-length(atom)-length(x) )
        atom=string( atom, x)
        # Y
        y=string( round(atoms.positions[i,2], digits=3 ) )
        atom=utils.spaces( atom, 46-length(atom)-length(y) )
        atom=string( atom, y )
        # Z
        z=string( round( atoms.positions[i,3], digits=3 ))
        atom=utils.spaces( atom, 54-length(atom)-length(z) )
        atom=string( atom, z )

        # Occupancy and Temperature Factor (default 0, except for PLUMED where used for masses)
        # Occupancy
        occ = periodicTable.names2Z( atoms.names[i] )
        atom=utils.spaces(atom, 60-length(atom)-length(occ) )
        atom=string(atom, occ )
        # Temperature Factor
        tempfac = periodicTable.names2Z( atoms.names[i] )
        atom=utils.spaces(atom,66-length(atom))
        atom=string(atom, string( 0.0 ),string(0) )

        # Atom name (right justified)
        atom=utils.spaces(atom,78-length(atom)-length(atoms.names[i]))
        atom=string(atom,atoms.names[i])

        # Charge (2 col max, default is blank)
        charge=string(" ")
        atom=utils.spaces(atom,80-length(atom)-length(charge))
        atom=string(atom,charge)

        # End of line
        atom=string(atom,"\n")
        Base.write(handle_out,atom)
    end

    Base.write(handle_out,"END\n")

    close(handle_out)

    return
end
function writePdb( file::T1, traj::Vector{T2}, cell::T3 ) where { T1 <: AbstractString , T2 <: atom_mod.AtomList, T3 <: cell_mod.Cell_param }

    # Writes a PDB file according to standard format (2011)

    handle_out=open(file,"w")

    nb_step=size(traj)
    for step=1:nb_step

        if ! writeCRYST1( handle_out, cell.length, cell.angles )
            print("Error writting cell information!\n")
            return false
        end

        # Atomic Positions
        nb_atoms = size(traj[step].names)[1]
        for i=1:nb_atoms
            atom="ATOM"
            # Right justified
            atom=utils.spaces(atom,11-length(string(traj[step].index[i]))-length(atom))
            atom=string(atom,traj[step].index[i])
            # If atom name is 1 length, start on 13, otherwise 14
            # atom names are left justified here
            if length(traj[step].names[i])==1
                atom=utils.spaces(atom,13-length(atom))
            else
                atom=utils.spaces(atom,14-length(atom))
            end
            atom=string(atom,traj[step].names[i])

            # Alternate location indicator ( default blank )
            alt_loc=string(" ")
            atom=utils.spaces(atom,17-length(atom)-length(alt_loc))
            atom=string(atom," ")

            # Residue name (3 col max) righ justified?
            residue_name=string("   ")
            atom=utils.spaces(atom,20-length(atom)-length(residue_name))
            atom=string(atom,residue_name)

            # Chain Identifier
            # Default is "X"
            chain_id=string("X")
            atom=utils.spaces(atom,22-length(atom)-length(chain_id))
            atom=string(atom,"X")

            # Residue Sequence Nb - Molecule nb
            mol_nb=string("1")
            atom=utils.spaces(atom,26-length(atom)-length(mol_nb))
            atom=string(atom,mol_nb)

            # Code for insertion of residues
            code_insertion=string(" ")
            atom=utils.spaces(atom,27-length(atom)-length(code_insertion))
            atom=string(atom,code_insertion)

            # Positions are right justified
            # X
            x=string( round(traj[step].positions[i,1], digits=3 ) )
            atom=utils.spaces( atom, 38-length(atom)-length(x) )
            atom=string( atom, x)
            # Y
            y=string( round(traj[step].positions[i,2], digits=3 ) )
            atom=utils.spaces( atom, 46-length(atom)-length(y) )
            atom=string( atom, y )
            # Z
            z=string( round( traj[step].positions[i,3], digits=3 ))
            atom=utils.spaces( atom, 54-length(atom)-length(z) )
            atom=string( atom, z )

            # Occupancy and Temperature Factor (default 0, except for PLUMED where used for masses)
            # Occupancy
            occ = periodicTable.names2Z( traj[step].names[i] )
            atom=utils.spaces(atom, 60-length(atom)-length(occ) )
            atom=string(atom, occ )
            # Temperature Factor
            tempfac = periodicTable.names2Z( traj[step].names[i] )
            atom=utils.spaces(atom,66-length(atom))
            atom=string(atom, string( 0.0 ),string(0) )

            # Atom name (right justified)
            atom=utils.spaces(atom,78-length(atom)-length(traj[step].names[i]))
            atom=string(atom,traj[step].names[i])

            # Charge (2 col max, default is blank)
            charge=string(" ")
            atom=utils.spaces(atom,80-length(atom)-length(charge))
            atom=string(atom,charge)

            # End of line
            atom=string(atom,"\n")
            Base.write(handle_out,atom)
        end

        Base.write(handle_out,"END\n")
    end

    close(handle_out)

    return
end
function writePdb( file::T1, traj::Vector{T2}, cells::Vector{T3} ) where { T1 <: AbstractString , T2 <: atom_mod.AtomList, T3 <: cell_mod.Cell_param }

    # Writes a PDB file according to standard format (2011)

    handle_out=open(file,"w")

    nb_step=size(traj)
    for step=1:nb_step

        if ! writeCRYST1( handle_out, cells[step].length, cells[step].angles )
            print("Error writting cell information!\n")
            return false
        end

        # Atomic Positions
        nb_atoms = size(traj[step].names)[1]
        for i=1:nb_atoms
            atom="ATOM"
            # Right justified
            atom=utils.spaces(atom,11-length(string(traj[step].index[i]))-length(atom))
            atom=string(atom,traj[step].index[i])
            # If atom name is 1 length, start on 13, otherwise 14
            # atom names are left justified here
            if length(traj[step].names[i])==1
                atom=utils.spaces(atom,13-length(atom))
            else
                atom=utils.spaces(atom,14-length(atom))
            end
            atom=string(atom,traj[step].names[i])

            # Alternate location indicator ( default blank )
            alt_loc=string(" ")
            atom=utils.spaces(atom,17-length(atom)-length(alt_loc))
            atom=string(atom," ")

            # Residue name (3 col max) righ justified?
            residue_name=string("   ")
            atom=utils.spaces(atom,20-length(atom)-length(residue_name))
            atom=string(atom,residue_name)

            # Chain Identifier
            # Default is "X"
            chain_id=string("X")
            atom=utils.spaces(atom,22-length(atom)-length(chain_id))
            atom=string(atom,"X")

            # Residue Sequence Nb - Molecule nb
            mol_nb=string("1")
            atom=utils.spaces(atom,26-length(atom)-length(mol_nb))
            atom=string(atom,mol_nb)

            # Code for insertion of residues
            code_insertion=string(" ")
            atom=utils.spaces(atom,27-length(atom)-length(code_insertion))
            atom=string(atom,code_insertion)

            # Positions are right justified
            # X
            x=string( round(traj[step].positions[i,1], digits=3 ) )
            atom=utils.spaces( atom, 38-length(atom)-length(x) )
            atom=string( atom, x)
            # Y
            y=string( round(traj[step].positions[i,2], digits=3 ) )
            atom=utils.spaces( atom, 46-length(atom)-length(y) )
            atom=string( atom, y )
            # Z
            z=string( round( traj[step].positions[i,3], digits=3 ))
            atom=utils.spaces( atom, 54-length(atom)-length(z) )
            atom=string( atom, z )

            # Occupancy and Temperature Factor (default 0, except for PLUMED where used for masses)
            # Occupancy
            occ = periodicTable.names2Z( traj[step].names[i] )
            atom=utils.spaces(atom, 60-length(atom)-length(occ) )
            atom=string(atom, occ )
            # Temperature Factor
            tempfac = periodicTable.names2Z( traj[step].names[i] )
            atom=utils.spaces(atom,66-length(atom))
            atom=string(atom, string( 0.0 ),string(0) )

            # Atom name (right justified)
            atom=utils.spaces(atom,78-length(atom)-length(traj[step].names[i]))
            atom=string(atom,traj[step].names[i])

            # Charge (2 col max, default is blank)
            charge=string(" ")
            atom=utils.spaces(atom,80-length(atom)-length(charge))
            atom=string(atom,charge)

            # End of line
            atom=string(atom,"\n")
            Base.write(handle_out,atom)
        end

        Base.write(handle_out,"END\n")
    end

    close(handle_out)

    return
end
#------------------------------------------------------------------------------#
function writePdbPlumed(atoms::T1, cell::T2, file::T3 ) where { T1 <: atom_mod.AtomMolList, T2 <: cell_mod.Cell_param, T3 <: AbstractString }

  out=open(file,"w")

  a,b,c = string(cell.length[1]), string(cell.length[2]), string(cell.length[3])
  alpha, beta, gamma = string(cell.angles[1]), string(cell.angles[2]), string(cell.angles[3])

  cryst1=string("CRYST1 ",a)
  cryst1=utils.spaces(cryst1,16-length(cryst1))
  cryst1=string(cryst1,b)
  cryst1=utils.spaces(cryst1,25-length(cryst1))
  cryst1=string(cryst1,c)
  cryst1=utils.spaces(cryst1,34-length(cryst1))
  cryst1=string(cryst1,alpha)
  cryst1=utils.spaces(cryst1,41-length(cryst1))
  cryst1=string(cryst1,beta)
  cryst1=utils.spaces(cryst1,48-length(cryst1))
  cryst1=string(cryst1,gamma)
  cryst1=utils.spaces(cryst1,56-length(cryst1))
  cryst1=string(cryst1,"P 1")
  cryst1=utils.spaces(cryst1,67-length(cryst1))
  cryst1=string(cryst1,"1")
  cryst1=string(cryst1,"\n")
  Base.write(out,cryst1)

  nb_atoms = size(atoms.atom_names)[1]
  for i=1:nb_atoms
    atom="ATOM"
    atom=utils.spaces(atom,11-length(string(atoms.atom_index[i]))-length(atom))
    atom=string(atom,atoms.atom_index[i])
    atom=utils.spaces(atom,13-length(atom))
    atom=string(atom,atoms.atom_names[i])
    atom=utils.spaces(atom,17-length(atom))
    atom=string(atom,atoms.mol_names[i])
    atom=utils.spaces(atom,21-length(atom))
    atom=string(atom,"X")
    atom=utils.spaces(atom,25-length(atom))
    atom=string(atom,atoms.mol_index[i])
    atom=utils.spaces(atom,32-length(atom))
    atom=string(atom,round(atoms.positions[i,1],digits=3 ) )
    atom=utils.spaces(atom,40-length(atom))
    atom=string(atom,round(atoms.positions[i,2],digits=3 ) )
    atom=utils.spaces(atom,48-length(atom))
    atom=string(atom,round(atoms.positions[i,3],digits=3 ) )
    atom=utils.spaces(atom,56-length(atom))
    atom=string(atom, string( round( periodicTable.names2Z(atoms.names[i]) ) ) )
    atom=utils.spaces(atom,62-length(atom))
    atom=string(atom, string( round( periodicTable.names2Z(atoms.names[i]) ) ) )
    atom=utils.spaces(atom,77-length(atom))
    atom=string(atom,"\n")
    Base.write(out,atom)
  end

  Base.write(out,"END\n")

  close(out)

  return
end
function writePdbPlumed(atoms::T1, cell::T2, file::T3 ) where { T1 <: atom_mod.AtomList, T2 <: cell_mod.Cell_param, T3 <: AbstractString }

  out=open(file,"w")

  a,b,c = string(round(cell.length[1],digits=2)), string(round(cell.length[2],digits=2)), string(round(cell.length[3],digits=2))
  alpha, beta, gamma = string(round(cell.angles[1],digits=2)), string(round(cell.angles[2],digits=2)), string(round(cell.angles[3],digits=2))

  cryst1=string("CRYST1 ",a)
  cryst1=utils.spaces(cryst1,16-length(cryst1))
  cryst1=string(cryst1,b)
  cryst1=utils.spaces(cryst1,25-length(cryst1))
  cryst1=string(cryst1,c)
  cryst1=utils.spaces(cryst1,34-length(cryst1))
  cryst1=string(cryst1,alpha)
  cryst1=utils.spaces(cryst1,41-length(cryst1))
  cryst1=string(cryst1,beta )
  cryst1=utils.spaces(cryst1,48-length(cryst1))
  cryst1=string(cryst1,gamma )
  cryst1=utils.spaces(cryst1,56-length(cryst1))
  cryst1=string(cryst1,"P 1")
  cryst1=utils.spaces(cryst1,67-length(cryst1))
  cryst1=string(cryst1,"1")
  cryst1=string(cryst1,"\n")
  Base.write(out,cryst1)

  nb_atoms = size(atoms.names)[1]
  for i=1:nb_atoms
    atom="ATOM"
    atom=utils.spaces(atom,11-length(string(atoms.index[i]))-length(atom))
    atom=string(atom,atoms.index[i])
    atom=utils.spaces(atom,13-length(atom))
    atom=string(atom,atoms.names[i])
    atom=utils.spaces(atom,17-length(atom))
    atom=string(atom," ")
    atom=utils.spaces(atom,21-length(atom))
    atom=string(atom,"X")
    atom=utils.spaces(atom,25-length(atom))
    atom=string(atom,1)
    atom=utils.spaces(atom,32-length(atom))
    atom=string(atom,round(atoms.positions[i,1], digits=3 ) )
    atom=utils.spaces(atom,40-length(atom))
    atom=string(atom,round(atoms.positions[i,2], digits=3 ) )
    atom=utils.spaces(atom,48-length(atom))
    atom=string(atom,round(atoms.positions[i,3], digits=3 ) )
    atom=utils.spaces(atom,56-length(atom))
    atom=string(atom, string( 0.0 ),string(0) )
    atom=utils.spaces(atom,62-length(atom))
    atom=string(atom, string( 0.0 ),string(0) )
    atom=utils.spaces(atom,76-length(atom))
    atom=string(atom,atoms.names[i])
    atom=string(atom,"\n")
    Base.write(out,atom)
  end

  Base.write(out,"END\n")

  close(out)

  return
end
#==============================================================================#

end
