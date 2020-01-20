module atom_mod

export Atom, AtomList, AtomMolList, switchAtoms, move

using utils

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

function getPositions( traj::Vector{T1} ) where { T1 <: AtomList }
    nb_step=size(traj)[1]
    nb_atoms=size(traj[1].names)[1]
    positions=zeros(nb_step,nb_atoms,3)
    for step=1:nb_step
        positions[step,:,:] = traj[step].positions[:,:]
    end
    return positions
end

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

end
