module pdb

export getNbSteps, readStep, writePDB, writePDBplumed

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
function readStep( file_path::T1 ) where { T1 <: AbstractString }
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
function writePDB( atoms::T1, cell::T2, file::T3 ) where { T1 <: atom_mod.AtomMolList, T2 <: cell_mod.Cell_param, T3 <: AbstractString }

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

  Base.write(out,string("MODEL X\n"))

  nb_atoms = size(atoms.atom_names)[1]
  for i=1:nb_atoms
    atom="ATOM"
    atom=utils.spaces(atom,7-length(atom))
    atom=string(atom,atoms.atom_index[i])
    atom=utils.spaces(atom,13-length(atom))
    atom=string(atom,atoms.atom_names[i])
    atom=utils.spaces(atom,23-length(atom))
    atom=string(atom,atoms.mol_names[i])
    atom=utils.spaces(atom,27-length(atom))
    atom=string(atom,atoms.mol_index[i])
    atom=utils.spaces(atom,31-length(atom))
    atom=string(atom,round(atoms.positions[i,1]*1000)/1000)
    atom=utils.spaces(atom,39-length(atom))
    atom=string(atom,round(atoms.positions[i,2]*1000)/1000)
    atom=utils.spaces(atom,47-length(atom))
    atom=string(atom,round(atoms.positions[i,3]*1000)/1000)
    atom=utils.spaces(atom,55-length(atom))
    atom=string(atom,"0.00")
    atom=utils.spaces(atom,61-length(atom))
    atom=string(atom,"0.00")
    atom=utils.spaces(atom,77-length(atom))
    atom=string(atom,atoms.atom_names[i])
    atom=string(atom,"\n")
    Base.write(out,atom)
  end

  Base.write(out,"END\n")

  close(out)

  return
end
function writePDB( atoms::T1, cell::T2, file::T3 ) where { T1 <: atom_mod.AtomList, T2 <: cell_mod.Cell_param, T3 <: AbstractString }

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

  Base.write(out,string("MODEL X\n"))

  nb_atoms = size(atoms.names)[1]
  for i=1:nb_atoms
    atom="ATOM"
    atom=utils.spaces(atom,7-length(atom))
    atom=string(atom,atoms.index[i])
    atom=utils.spaces(atom,13-length(atom))
    atom=string(atom,atoms.names[i])
    atom=utils.spaces(atom,23-length(atom))
    atom=string(atom,"XXX")
    atom=utils.spaces(atom,27-length(atom))
    atom=string(atom,"XXX")
    atom=utils.spaces(atom,31-length(atom))
    atom=string(atom,round(atoms.positions[i,1]*1000)/1000)
    atom=utils.spaces(atom,39-length(atom))
    atom=string(atom,round(atoms.positions[i,2]*1000)/1000)
    atom=utils.spaces(atom,47-length(atom))
    atom=string(atom,round(atoms.positions[i,3]*1000)/1000)
    atom=utils.spaces(atom,55-length(atom))
    atom=string(atom,"0.00")
    atom=utils.spaces(atom,61-length(atom))
    atom=string(atom,"0.00")
    atom=utils.spaces(atom,77-length(atom))
    atom=string(atom,atoms.names[i])
    atom=string(atom,"\n")
    Base.write(out,atom)
  end

  Base.write(out,"END\n")

  close(out)

  return
end
function writePDB( atoms::T1, cell::T2, fileio::T3 ) where { T1 <: atom_mod.AtomList, T2 <: cell_mod.Cell_param, T3 <: IO }

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
  Base.write(fileio,cryst1)

  Base.write(fileio,string("MODEL X\n"))

  nb_atoms = size(atoms.names)[1]
  for i=1:nb_atoms
    atom="ATOM"
    atom=utils.spaces(atom,7-length(atom))
    atom=string(atom,atoms.index[i])
    atom=utils.spaces(atom,13-length(atom))
    atom=string(atom,atoms.names[i])
    atom=utils.spaces(atom,23-length(atom))
    atom=string(atom,"XXX")
    atom=utils.spaces(atom,27-length(atom))
    atom=string(atom,"XXX")
    atom=utils.spaces(atom,31-length(atom))
    atom=string(atom,round(atoms.positions[i,1]*1000)/1000)
    atom=utils.spaces(atom,39-length(atom))
    atom=string(atom,round(atoms.positions[i,2]*1000)/1000)
    atom=utils.spaces(atom,47-length(atom))
    atom=string(atom,round(atoms.positions[i,3]*1000)/1000)
    atom=utils.spaces(atom,55-length(atom))
    atom=string(atom,"0.00")
    atom=utils.spaces(atom,61-length(atom))
    atom=string(atom,"0.00")
    atom=utils.spaces(atom,77-length(atom))
    atom=string(atom,atoms.names[i])
    atom=string(atom,"\n")
    Base.write(fileio,atom)
  end

  Base.write(fileio,"END\n")

  return
end
function writePDBplumed(atoms::T1, cell::T2, file::T3 ) where { T1 <: atom_mod.AtomMolList, T2 <: cell_mod.Cell_param, T3 <: AbstractString }

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
function writePDBplumed(atoms::T1, cell::T2, file::T3 ) where { T1 <: atom_mod.AtomList, T2 <: cell_mod.Cell_param, T3 <: AbstractString }

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

function writeCRYST1( handle_out::T1, lenghts::Vector{T1}, angles::Vector{T2} ) where { T1 <: IO, T2 <: Real, T3 <: Real }
    return true
end

function writePdbVmd( file::T1, atoms::T2, cell::T3 ) where { T1 <: AbstractString , T2 <: atom_mod.AtomList, T3 <: cell_mod.Cell_param }

  handle_out=open(file,"w")

  a,b,c = string(round(cell.length[1],digits=3)), string(round(cell.length[2],digits=3)), string(round(cell.length[3],digits=3))
  alpha, beta, gamma = string(round(cell.angles[1],digits=2)), string(round(cell.angles[2],digits=2)), string(round(cell.angles[3],digits=2))
  space_group=string("P 1") # In Hermann Mauguin notation, by default P 1 works... never seen a case where something else is used 28/01/2020 (M Moog)
  z_number = string("1") # Number of polumeric chains or in hetero polymers the n umber of occurences of the most populous chain

  #Implement if longer cut?

  # Writting CRYST1 info -> Cell information
  cryst1=string("CRYST1")
  cryst1=utils.spaces(cryst1,15-length(cryst1)-length(a))
  cryst1=string(cryst1,a)
  cryst1=utils.spaces(cryst1,24-length(cryst1)-length(b))
  cryst1=string(cryst1,b)
  cryst1=utils.spaces(cryst1,33-length(cryst1)-length(c))
  cryst1=string(cryst1,c)
  cryst1=utils.spaces(cryst1,40-length(cryst1)-length(alpha))
  cryst1=string(cryst1,alpha)
  cryst1=utils.spaces(cryst1,47-length(cryst1)-length(beta))
  cryst1=string(cryst1,beta)
  cryst1=utils.spaces(cryst1,54-length(cryst1)-length(gamma))
  cryst1=string(cryst1,gamma )
  cryst1=utils.spaces(cryst1,56-length(cryst1))
  cryst1=string(cryst1,space_group)
  cryst1=utils.spaces(cryst1,70-length(cryst1)-length(z_number))
  cryst1=string(cryst1,z_number)
  cryst1=string(cryst1,"\n")
  Base.write(handle_out,cryst1)

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
#==============================================================================#

end
