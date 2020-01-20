module pdb

export getNbSteps, readStep, writePDB, writePDBplumed

using utils
using atom_mod
using cell_mod

function getNbSteps( file::T1 ) where { T1 <: AbstractString}
  count=0
  count2=0
  open( file ) do f
    while !eof(f)
      col1 =split(readline(f))[1]
      print("col1: ", col1, "\n")
      if col1 == "END"
        count += 1
      elseif col1 == "CRYST1"
        count2 +=1
      end
    end
  end
  print(count2)
  if count == count2
    return count
  else
    print("Problem in file...\n")
    return 0
  end
end

function readStep( file::T1 ) where { T1 <: AbstractString }
  #--------------
  # Reading file
  #----------------------
  file=open(file);
  lines=readlines(file);
  close(file);
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

#-------------------------------------------------------------------------------

function writePDB(atoms::T1, cell::T2, file::T3 ) where { T1 <: atom_mod.AtomMolList, T2 <: cell_mod.Cell_param, T3 <: AbstractString }

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
    atom=string(atom,round(atoms.positions[i,1]*1000)/1000)
    atom=utils.spaces(atom,40-length(atom))
    atom=string(atom,round(atoms.positions[i,2]*1000)/1000)
    atom=utils.spaces(atom,48-length(atom))
    atom=string(atom,round(atoms.positions[i,3]*1000)/1000)
    atom=utils.spaces(atom,56-length(atom))
    atom=string(atom,"8.00")
    atom=utils.spaces(atom,62-length(atom))
    atom=string(atom,"8.00")
    atom=utils.spaces(atom,77-length(atom))
    atom=string(atom,"\n")
    Base.write(out,atom)
  end

  Base.write(out,"END\n")

  close(out)

  return
end


function writePDB(atoms::T1, cell::T2, file::T3 ) where { T1 <: atom_mod.AtomList, T2 <: cell_mod.Cell_param, T3 <: AbstractString }

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

end
