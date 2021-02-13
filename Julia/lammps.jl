module cpmd

# Loading necessary modules from LibAtomSim
using conversion
using atom_mod

# Description
# Handles interface with LAMMPS: Reads input and output files, writes input file,
# sets up simulations based on user parameters, etc...

export writeLammpsInput

# Writting from a collection of molecules
#-------------------------------------------------------------------------------
function writeLammpsInput( molecules::Vector{T1}, file::T2 ) where { T1 <: atom_mod.AtomMolList, T2 <: AbstractString }
  # Variables
  # - number of atomic species
  nb_types=size(getNames(molecules))[1]
  # - number of atoms
  nb_atoms=getAtomsNb(molecules)
  # Opening File
  out=open(file,"w+")
  # Introduction
  #--------------------------------------------------
  write(out,"LAMMPS data file, timestep = 0\n")
  write(out,"\n")
  #--------------------------------------------------
  # Number of atoms and types
  #--------------------------------------
  write(out,string(nb_atoms, " atoms\n"))
  write(out,string(nb_types, " atom types\n"))
  #---------------------------------------
  # Delimiting the cell
  #-----------------------------------------------------------
  write(out,"\n")
  X=getX(molecules)
  write(out,string(getMin(X)," ",getMax(X)," xlo xhi\n"))
  empty!(X)
  Y=getY(molecules)
  write(out,string(getMin(Y)," ",getMax(Y)," ylo yhi\n"))
  empty!(Y)
  Z=getZ(molecules)
  write(out,string(getMin(Z)," ",getMax(Z)," zlo zhi\n"))
  empty!(Z)
  #-----------------------------------------------------------
  # Masses
  #---------------------------------------------
  write(out,"\n")
  write(out,"Masses\n")
  write(out,"\n")
  masses=getMasses(molecules)
  labels=getLabels(molecules)
  for i=1:nb_types
    write(out,string(labels[i]," ",masses[i],"\n"))
  end
  #----------------------------------------------
  write(out,"\n")
  write(out,"Atoms\n")
  write(out,"\n")
  # style of atoms in data file: atom-ID molecule-ID atom-type q x y z
  atom_nb=1
  for i=1:size(molecules)[1]
    atoms = getAtoms(molecules[i])
    for j=1:size(atoms)[1]
      write(out,string(atom_nb," "))
      write(out,string(i," "))
      write(out,string(getLabel(atoms[j])," "))
      write(out,string(getCharge(atoms[j])," "))
      write(out,string(getX(atoms[j])," ",getY(atoms[j])," ",getZ(atoms[j]) ))
      write(out,"\n")
      atom_nb=atom_nb+1
    end
  end
  close(out)
end
#-------------------------------------------------------------------------------

end
