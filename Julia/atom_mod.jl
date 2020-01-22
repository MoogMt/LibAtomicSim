module atom_mod

export Atom, AtomList, AtomMolList,
export switchAtoms, moveAtom, moveAtom!
export getNbAtoms, getNbStep, getNbMol

using utils

#-------------------------------------------------------------------------------
mutable struct Atom
    name::AbstractString
    index::Int
    position::Vector{Real}
    function Atom()
        new("",-1,[0,0,0]);
    end
end
mutable struct AtomList
    names::Vector{AbstractString}
    index::Vector{Int}
    positions::Array{Real}
    function AtomList( nb_atoms::T1 ) where { T1 <: Int }
        new( Array{AbstractString,1}(undef,nb_atoms) , zeros(nb_atoms), zeros(nb_atoms,3) )
    end
end
mutable struct AtomMolList
    atom_names::Vector{AbstractString}
    atom_index::Vector{Int}
    mol_names::Vector{AbstractString}
    mol_index::Vector{Int}
    positions::Array{Real}
    function AtomMolList( nb_atoms::T1 ) where { T1 <: Int }
        new( Vector{AbstractString}( undef,nb_atoms ),Vector{Int}(undef, nb_atoms ) ,Vector{AbstractString}(undef, nb_atoms ), Vector{Int}(undef,nb_atoms), Array{Real}(undef,nb_atoms,3))
    end
end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
function switchAtoms!( atoms::T1 , index1::T2, index2::T3 ) where { T1 <: AtomList, T2 <: Int, T3 <: Int }
    # Storing
    index_    = atoms.index[ index1 ]
    name_     = atoms.names[ index1 ]
    position_ = atoms.positions[ index1, : ]
    # Moving 1
    atoms.index[ index1 ]  = atoms.index[ index2 ]
    atoms.names[ index1 ]  = atoms.names[ index2 ]
    atoms.positions[ index1, : ] = atoms.positions[ index2, : ]
    # Moving 2
    atoms.index[index2]       = index_
    atoms.names[index2]       = name_
    atoms.positions[index2,:] = position_
    return
end
function switchAtoms!( traj::Vector{T1} , index1::T2, index2::T3, step::T4 ) where { T1 <: AtomList, T2 <: Int, T3 <: Int, T4 <: Int }
    return SwitchAtoms!( traj[step], index1, index2 )
end
function switchAtoms!( atoms::T1 , index1::T2, index2::T3 ) where { T1 <: AtomMolList, T2 <: Int, T3 <: Int }
    # Storing
    a_index=atoms.atom_index[index1]
    a_name=atoms.atom_names[index1]
    m_index=atoms.mol_index[index1]
    m_name=atoms.mol_names[index1]
    positions=atoms.positions[index1,:]
    # Moving 1
    atoms.atom_index[index1]=atoms.atom_index[index2]
    atoms.atom_names[index1]=atoms.atom_names[index2]
    atoms.mol_index[index1]=atoms.mol_index[index2]
    atoms.mol_names[index1]=atoms.mol_names[index2]
    atoms.positions[index1,:]=atoms.positions[index2,:]
    # Moving 2
    atoms.atom_index[index2]=a_index
    atoms.atom_names[index2]=a_name
    atoms.mol_index[index2]=m_index
    atoms.mol_names[index2]=m_name
    atoms.positions[index2,:]=positions
    return
end
function switchAtoms!( traj::Vector{T1}, index1::T2, index2::T3, step::T4 ) where { T1 <: AtomMolList, T2 <: Int, T3 <: Int, T4 <: Int }
    return switchAtoms!( traj[step], index1, index2 )
end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
function moveAtom( atoms::T1 , move::Vector{T2} ) where { T1 <: atom_mod.AtomList, T2 <: Real }
    for i=1:size(atoms.positions)[1]
        for j=1:3
            atoms.positions[i,j] += move[j]
        end
    end
    return atoms
end
function moveAtom!( atoms::T1 , move::Vector{T2} ) where { T1 <: atom_mod.AtomList, T2 <: Real }
    for i=1:size(atoms.positions)[1]
        for j=1:3
            atoms.positions[i,j] += move[j]
        end
    end
    return atoms
end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
function getPositionsAsArray( traj::Vector{T1} ) where { T1 <: AtomList }
    nb_step=size(traj)[1]
    nb_atoms=size(traj[1].names)[1]
    positions=zeros(nb_step,nb_atoms,3)
    for step=1:nb_step
        positions[step,:,:] = traj[step].positions[:,:]
    end
    return positions
end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
function getNbAtoms( atoms::T1 ) where { T1 <: AtomList}
    return size( atoms.names )[1]
end
function getNbAtoms( traj::Vector{T1} ) where { T1 <: AtomList }
    return getNbAtoms( traj[1] )
end
function getNbAtoms( atoms::T1 ) where { T1 <: AtomMolList }
    return size( atoms.atom_names )[1]
end
function getNbAtoms( traj::Vector{T1} ) where { T1 <: AtomMolList}
    return getNbAtoms( traj[1] )
end
function getNbMolecules( atoms::T1 ) where { T1 <: AtomMolList }
    return max( atoms.mol_index )
end
function getNbMolecules( traj::Vector{T1} ) where { T1 <: AtomMolList }
    return getNbMolecules( traj[1] )
end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
function getNbStep( traj::Vector{T1} ) where { T1 <: AtomList }
    return size(traj)[1]
end
function getNbStep( traj::Vector{T1} ) where { T1 <: AtomMolList }
    return size(traj)[1]
end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
function getTypeIndex( types_names::Vector{T1}, name::T2 ) where { T1 <: AbstractString , T2 <: AbstractString }
    index_types=zeros(Int,0)
    nb_atoms=size(types_names)[1]
    for atom=1:nb_atoms
        if types_names[atom] == name
            push!(index_types,atom)
        end
    end
    return index_types
end
function getSpecies( atoms::T1 ) where { T1 <: AtomList }
    species = Vector{AbstractString}(undef,0)
    nb_atoms = getNbAtoms( atoms )
    if nb_atoms == 0
        return []
    else
        push!( species, atoms.names[1] )
        for atom = 1:nb_atoms
            add = true
            for specie in species
                if atoms.names[atom] == specie
                    add = false
                    break
                end
            end
            if add
                push!( species, atom.names[atom] )
            end
        end
        return species
    end
end
function getSpecies( traj::Vector{T1} ) where { T1 <: AtomList }
    nb_step = getNbStep( traj )
    species=[]
    for step = 1:nb_step
        push!(species, getSpecies( traj[step] ) )
    end
    return species
end
function getNbSpecies( atoms::T1 ) where { T1 <: AtomList }
    nb_atoms = getNbAtoms( atoms )
    if nb_atoms == 0
        return 0
    else
        return size(getSpecies( atoms ))[1]
    end
end
function getStartSpecie( atoms::T1, species::Vector{T2} ) where { T1 <: AtomList, T2 <: AbstractString }
#-------------------------------------------------------------------------------

end
