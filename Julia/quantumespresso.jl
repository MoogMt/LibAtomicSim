module quantumespresso

#===================#
# QE Output reading #
#===================#

#--------#
# Energy
#------------------------------------------------------------------------------#
# Get the energy from a simulation
function getEnergyQE(fname::AbstractString)
  energy=[]
  open(fname,"r") do f
    for line in eachline(f)
      if( contains(line,"!") )
        push!(energy,parse(Float64,AbstractString(split(split(line,"=")[2],"R")[1])))
      end
    end
  end
  return energy
end
# Relax Energy
function relaxEnergy(fname::AbstractString)
  energy=getEnergyQE(fname)
  return energy[size(energy)[1]]
end
# get relax binding energy
function bindingE{T1 <: Int, T2 <: Real}(fname::AbstractString,specie_number::Vector{T1},isolated_energy::Vector{T2})
  size1=size(specie_number)[1]
  if size1 != size(isolated_energy)[1]
    error("Vector sizes do not match")
    return
  end
  atoms_nb=0
  energy=relaxEnergy(fname)
  for i=1:size1
    energy = energy - specie_number[i]*isolated_energy[i]
    atoms_nb = atoms_nb + specie_number[i]
  end
  return convertRy2eV(-energy/atoms_nb)
end
#------------------------------------------------------------------------------#

#-
# Magnetization
#------------------------------------------------------------------------------#
function getMagQE(file::AbstractString)
  mag=Array{Float64}(0)
  open(file,"r") do f
    for line in eachline(f)
      if contains(line,"total magnetization")
        push!(mag,parse(Float64,AbstractString(split(split(line,"=")[2],"Bohr")[1])))
      end
    end
  end
  return mag
end
function relaxMag(file::AbstractString)
  mag=getMagQE(file)
  return mag[size(mag)[1]]
end
#------------------------------------------------------------------------------#

end
