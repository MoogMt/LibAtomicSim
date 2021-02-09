module atom_mod

export AtomList, AtomMolList
export switchAtoms!, moveAtom! sortAtomsByZ!
export getNbAtoms, getNbStep, getNbMol, moveAtom

using utils
using periodicTable

# Structures
#-------------------------------------------------------------------------------
# AtomList: Contains all atoms with positions, names and index
mutable struct AtomList

    # Variables
    #-------------------------------
    names::Vector{AbstractString} # Names of the atoms (chemical species)
    index::Vector{Int}            # Index of atoms (labels)
    positions::Array{Real,2}      # Positions of the atoms
    #-------------------------------

    # Constructors
    #----------------------------------------------------------------------------
    # Creates a default AtomList
    function AtomList()
        # Arguments
        # None
        # Output
        # - Creates a default AtomList

        # Create the new AtomList with nothing in it
        new( Array{AbstractString,1}(undef,0) , zeros(0), zeros(0,3) )
    end
    # Creates a default AtomList with a preset number of atoms
    function AtomList( nb_atoms::T1 ) where { T1 <: Int }
        # Arguments:
        # - nb_atoms : number of atoms in the AtomList
        # Output:
        # Creates an AtomList with nb_atoms, all with default values

        # Creates the AtomList
        new( Array{AbstractString,1}(undef,nb_atoms) , zeros(nb_atoms), zeros(nb_atoms,3) )
    end
    # Creates an AtomList with all data given
    function AtomList( names::Vector{T1}, index::Vector{T2}, positions::Array{T3,2} ) where { T1 <: AbstractString, T2 <: Int, T3 <: Real }
        # Arguments
        # - names: names of the atoms (vector of string)
        # - index: indexes of the atoms (vector of int)
        # - positions: array with the positions of the atoms (nb_atoms,3)
        # Output
        # - Creates an AtomList

        # Create the AtomList
        new( names, index, positions )
    end
    #----------------------------------------------------------------------------
end
# AtomMolList: Contains all atoms with positions, names and index + name and index of the molecule they belong to
# Mostly useful for classical simulations (GROMACS,LAMMPS)
mutable struct AtomMolList

    # Variables
    #---------------------------------------------------
    atom_names::Vector{AbstractString} # Atoms names
    atom_index::Vector{Int}            # Atom indexes
    mol_names::Vector{AbstractString}  # Names of the molecule the atom is in
    mol_index::Vector{Int}             # Indexes of the molecule the atom belong to
    positions::Array{Real,2}           # Position of the atom
    #---------------------------------------------------

    # Constructors
    #---------------------------------------------------
    # Create default AtomMolList
    function AtomMolList()
        # Argument
        # - None
        # Output
        # - A default AtomMolList object

        # Creates object
        new( Vector{AbstractString}( undef, 0 ),Vector{Int}(undef, 0 ) ,Vector{AbstractString}(undef, 0 ), Vector{Int}(undef, nb_atoms ), Array{Real}(undef, 0, 3 ) )
    end
    # Create default AtomMolList with a given number of atoms
    function AtomMolList( nb_atoms::T1 ) where { T1 <: Int }
        # Argument
        # - nb_atoms: number of atoms of the AtomMolList object
        # Output
        # - An AtomMolList object
        new( Vector{AbstractString}(undef, nb_atoms ),Vector{Int}(undef, nb_atoms ) ,Vector{AbstractString}(undef, nb_atoms ), Vector{Int}(undef, nb_atoms ), Array{Real}(undef, nb_atoms, 3 ) )
    end
    # Create AtomMolList
    function AtomMolList( atom_names::Vector{T1}, atom_index::Vector{T2}, mol_names::Vector{T3}, mol_index::Vector{T4}, positions::Array{T5,2} ) where { T1 <: AbstractString, T2 <: Int, T3 <: AbstractString, T4 <: Int, T5 <: Real }
        # Argument
        # - atom_names: names of the atoms (vector of Strings)
        # - atom_index: indexes of the atoms (vector of int)
        # - mol_names: names of the molecule the atoms belong to (Vector of string)
        # - mol_index: indexes of the molecules the atoms are in (Vector of int)
        # - positions: position of the atom (Array, real, (nb_atoms,3) )
        # Output
        # - An AtomMolList object

        # Creates the AtomMolList object
        new( Vector{AbstractString}(undef, nb_atoms ),Vector{Int}(undef, nb_atoms ) ,Vector{AbstractString}(undef, nb_atoms ), Vector{Int}(undef, nb_atoms ), Array{Real}(undef, nb_atoms, 3 ) )
    end
    #---------------------------------------------------
end
#-------------------------------------------------------------------------------

# Functions to switch atoms
#-------------------------------------------------------------------------------
# Switch two atoms in an AtomList
function switchAtoms!( atoms::T1 , index1::T2, index2::T3 ) where { T1 <: AtomList, T2 <: Int, T3 <: Int }
    # Argument
    # - atoms: AtomList containing information about the atoms and their positions
    # - index1, index2: index of the atoms to switch (int)
    # Output
    # - None

    # Stores data from 1
    index_    = atoms.index[ index1 ]
    name_     = atoms.names[ index1 ]
    position_ = atoms.positions[ index1, : ]

    # Moving data from 2 to 1
    atoms.index[ index1 ]        = atoms.index[ index2 ]
    atoms.names[ index1 ]        = atoms.names[ index2 ]
    atoms.positions[ index1, : ] = atoms.positions[ index2, : ]

    # Moving data from storage to2 2
    atoms.index[index2]       = index_
    atoms.names[index2]       = name_
    atoms.positions[index2,:] = position_

    return
end
# Switch two atoms from a traj (vector of AtomList)
function switchAtoms!( traj::Vector{T1}, index1::T2, index2::T3, step::T4 ) where { T1 <: AtomList, T2 <: Int, T3 <: Int, T4 <: Int }
    # Argument
    # - traj: vector of AtomList with positions of atoms
    # - index1, index2: index of the atoms
    # - step: target set
    # Output
    # - None

    # Switches the atom index1 and index2
    return SwitchAtoms!( traj[step], index1, index2 )
end
# Switches to atoms in an AtomMolList
function switchAtoms!( atoms::T1 , index1::T2, index2::T3 ) where { T1 <: AtomMolList, T2 <: Int, T3 <: Int }
    # Argument
    # - atoms:
    # - index1, index2: index of atom to be switched
    # Output
    # - None

    # Stores data from atom index 1
    a_index   = atoms.atom_index[index1]
    a_name    = atoms.atom_names[index1]
    m_index   =  atoms.mol_index[index1]
    m_name    =  atoms.mol_names[index1]
    positions =  atoms.positions[index1,:]

    # Moves data from atom index2 into atom index1
    atoms.atom_index[index1]  = atoms.atom_index[index2]
    atoms.atom_names[index1]  = atoms.atom_names[index2]
    atoms.mol_index[index1]   = atoms.mol_index[index2]
    atoms.mol_names[index1]   = atoms.mol_names[index2]
    atoms.positions[index1,:] = atoms.positions[index2,:]

    # Moves data from storage to atom index2
    atoms.atom_index[index2]  = a_index
    atoms.atom_names[index2]  = a_name
    atoms.mol_index[index2]   = m_index
    atoms.mol_names[index2]   = m_name
    atoms.positions[index2,:] = positions

    return
end
# Switches two atoms from traj (vector of AtomMolList)
function switchAtoms!( traj::Vector{T1}, index1::T2, index2::T3, step::T4 ) where { T1 <: AtomMolList, T2 <: Int, T3 <: Int, T4 <: Int }
    # Argument
    # - traj: vector of AtomMolList, contains atom positions
    # - index1, index2 : index of the atoms to switchs
    # - step: step at which to exchange atoms
    # Output
    # None

    # Switches atoms index1 and index2 at step
    return switchAtoms!( traj[step], index1, index2 )
end
#-------------------------------------------------------------------------------

# Sorting atoms by chemical species
#-------------------------------------------------------------------------------
# Sort atoms by their chemical species
function sortAtomsByZ!( atoms::T1 ) where { T1 <: AtomList }
    # Argument
    # - atoms: AtomList containing atomic positions
    # Output
    # - None

    # Get the number of atoms
    nb_atoms = size(atoms.names)[1]

    # Loop over atoms
    for atom = 1:nb_atoms
        # Loop over all other atoms
        for atom2 = atom+1:nb_atoms
            # Get and compare the Z of the two atoms
            if periodicTable.names2Z( atoms.names[atom] ) < periodicTable.names2Z( atoms.names[atom2] )
                # Switch atoms
                switchAtoms!( atoms, atom, atom2 )
            end
        end
    end

    return
end
# Sort atoms by their chemical species over all steps
function sortAtomsByZ!( traj::Vector{T1} ) where { T1 <: AtomList }
    # Argument
    # - traj: trajectory as list of AtomList containing the positions
    # Output
    # - None

    # Get the number of step
    nb_step = size(traj)[1]

    # Loop over steps
    for step = 1:nb_step
        # Sort for each step
        sortAtomsByZ!( traj[step] )
    end

    return
end
#-------------------------------------------------------------------------------

# Moving atoms around
#-------------------------------------------------------------------------------
# Moving all atoms by a given vector move
function moveAllAtom( atoms::T1 , move::Vector{T2} ) where { T1 <: atom_mod.AtomList, T2 <: Real }
    # Argument
    # - atoms: AtomList, with atomic positions
    # - move: vector by which to moive all atoms
    # Output
    # - atoms2, transformed: AtomList with the modified positions

    # Initialize output
    atoms2 = copy(atoms)

    # Loop over atoms
    for atom=1:size(atoms.positions)[1]
        # Loop over dimensions
        for i=1:3
            atoms2.positions[atom,i] += move[i]
        end
    end

    # Returns the modified atomic positions
    return atoms2
end
# Moving all atoms by a given vector, directly on the vector
function moveAllAtom!( atoms::T1 , move::Vector{T2} ) where { T1 <: atom_mod.AtomList, T2 <: Real }
    # Argument
    # - atoms: frame with atomic position described by AtomList
    # - move: vector by which moving all atoms of the list
    # Output
    # None (atoms are directly modified)

    # Loop over atoms
    for atom=1:size(atoms.positions)[1]
        # Loop over dimensions
        for i=1:3
            atoms.positions[atom,i] += move[i]
        end
    end

    return
end
#-------------------------------------------------------------------------------

# Extract positions of an AtomList
#-------------------------------------------------------------------------------
function getPositionsAsArray( traj::Vector{T1} ) where { T1 <: AtomList }
    # Argument
    # - traj: vector of AtomList that describes the trajectory
    # Output
    # - positions: Array containing atomic positions real, size (nb_atoms,3)

    # Get number of steps of the trajectory
    nb_step=size(traj)[1]

    # Get the number of atoms in the structure
    nb_atoms=size(traj[1].names)[1]

    # Initialize positions
    positions=zeros(nb_step,nb_atoms,3)

    # Loop over step
    for step=1:nb_step
        # Extracting position at each step
        positions[step,:,:] = traj[step].positions[:,:]
    end

    # Returns positions
    return positions
end
#-------------------------------------------------------------------------------

# Functions to get the number of steps/atoms from
# - an AtomList or vector of AtomList
# - an AtomMolList or vector of AtomMolList
#-------------------------------------------------------------------------------
# Extract the number of atoms from an AtomList
function getNbAtoms( atoms::T1 ) where { T1 <: AtomList }
    # Argument
    # - atoms: Structure described by an AtomList
    # Return
    # - number of atoms in the AtomList

    # Returns the number of atoms in AtomList
    return size( atoms.names )[1]
end
# Extract the number of atoms from a vector of AtomList
# - assumes that this is a trajectory, not a set of structure with different number of atoms
function getNbAtoms( traj::Vector{T1} ) where { T1 <: AtomList }
    # Argument
    # - traj: vector of AtomList, containing atom positions
    # Output
    # - number of atoms in the trajectory

    # Return the number of atoms in the trajectory
    return getNbAtoms( traj[1] )
end
# Extract the number of atoms from an AtomMolList
function getNbAtoms( atoms::T1 ) where { T1 <: AtomMolList }
    # Argument
    # - atoms: Structure described by an AtomMolList
    # Return
    # - number of atoms in the AtomMolList

    # Returns the number of atoms in AtomMolList
    return size( atoms.atom_names )[1]
end
# Extract the number of atoms from a vector of AtomMolList
# - assumes that this is a trajectory, not a set of structure with different number of atoms
function getNbAtoms( traj::Vector{T1} ) where { T1 <: AtomMolList }
    # Argument
    # - traj: vector of AtomMolList, containing atom positions
    # Output
    # - number of atoms in the trajectory

    # Return the number of atoms in the trajectory
    return getNbAtoms( traj[1] )
end
# Get the maximum number of molecules in an AtomMolList
function getNbMolecules( atoms::T1 ) where { T1 <: AtomMolList }
    # Argument
    # - atoms: AtomMolList containing atomic positions
    # Output
    # - number of molecules in the frame

    # Return the number of molecules (maximum index of molecules)
    return max( atoms.mol_index )
end
# Get the maximum number of molecules in a vector of AtomMolList
function getNbMolecules( traj::Vector{T1} ) where { T1 <: AtomMolList }
    # Argument
    # - traj: vector of AtomMolList containing atomic positions of the trajectory
    # Output
    # - number of molecules in the trajectory (assuming it does not change)

    # Return the number of molecules (maximum index of molecules)
    return getNbMolecules( traj[1] )
end
# Get the number of steps in the trajectory made up of a vector of AtomList
function getNbStep( traj::Vector{T1} ) where { T1 <: AtomList }
    # Argument
    # - traj: vector of AtomList
    # Output
    # number of steps of the vector

    # Return the size of the trajectory
    return size( traj )[1]
end
# Get the number of steps in the trajectory made up of a vector of AtomMolList
function getNbStep( traj::Vector{T1} ) where { T1 <: AtomMolList }
    # Argument
    # - traj: vector of AtomMolList
    # Output
    # number of steps of the vector

    # Return the size of the trajectory
    return size( traj )[1]
end
#-------------------------------------------------------------------------------

# Get information about the chemical species and type
#-------------------------------------------------------------------------------
# Get all indexes for a given specie
function getTypeIndex( types_names::Vector{T1}, name::T2 ) where { T1 <: AbstractString , T2 <: AbstractString }
    # Argument
    # - type_names: types of all the atoms in the structure
    # - name: type of atoms one wants the indexes of
    # Output
    # - index_types: list of index of atoms of type "name"

    # Initialize output
    index_types=zeros(Int,0)

    # Get the number of atoms in the list
    nb_atoms=size(types_names)[1]

    # Loop over atoms
    for atom=1:nb_atoms
        # If the atom if of the right type, add the index to the list
        if types_names[atom] == name
            push!( index_types, atom )
        end
    end

    # Get all indexes of the atoms
    return index_types
end
# Get all the species names within the AtomList frame
function getSpecies( atoms::T1 ) where { T1 <: AtomList }
    # Argument
    # - atoms: AtomList describing the structure
    # Output
    # - a list of all the species present in the AtomList

    # Initialize output
    species = Vector{AbstractString}(undef,0)

    # Get number of atoms
    nb_atoms = getNbAtoms( atoms )

    # If the number atom is null...
    if nb_atoms == 0
        # Returns empty list
        return []
    else
        # Add the type of the first atom to the list
        push!( species, atoms.names[1] )

        # Loop over the atoms
        for atom=2:nb_atoms

            # Check that the type already exists
            add = true

            # Loop over the species
            for specie in species
                # If type already exists in the list
                if atoms.names[atom] == specie
                    # Break and moves to the next atom
                    add = false
                    break
                end
            end

            # If the type was not found, adding it to the list
            if add
                push!( species, atoms.names[atom] )
            end
        end

        # Return the list of chemical species
        return species
    end
end
# Get all the species names within a vector of AtomList frame
# NB delete redundancy
function getSpecies( traj::Vector{T1} ) where { T1 <: AtomList }
    # Argument
    # - traj: vector of AtomList
    # Output
    # - vector of the names of all types of atom in the traj

    # Get the number of steps of the traj
    nb_step = getNbStep( traj )

    # Initialize output
    species = Vector{AbstractString}(undef,0)

    # Loop over steps
    for step = 1:nb_step
        # add all the species from each step (this means there will be redundancy...)
        push!( species, getSpecies( traj[step] ) )
    end

    # Return a vector with all species names in the traj
    return species
end
# Get the number of species in the AtomList
function getNbSpecies( atoms::T1 ) where { T1 <: AtomList }
    # Argument
    # - atoms: AtomList frame describing the structure
    # Output
    # - The number of species (int)

    # Get the number of atoms of AtomList
    nb_atoms = getNbAtoms( atoms )

    # If no atoms, there is 0 species
    if nb_atoms == 0
        return 0
    else
        # Compute and get the number of species
        return size( getSpecies( atoms ) )[1]
    end
end
# Get the number of species in vector of AtomList
function getNbSpecies( traj::Vector{T1} ) where { T1 <: AtomList }
    # Argument
    # - traj: vector of AtomList that contains the trajectory (atomic positions, etc...)
    # Output
    # - Number of species for each step (int)

    # Get the number of step in the trajectory
    nb_step = getNbStep( traj )

    # Initialize output
    species_nb = zeros( Int, nb_step )

    # Loop over steps
    for step = 1:nb_step
        # Get number of species for each step
        species_nb[ step ] = getNbSpecies( traj[step] )
    end

    # Return vector of number of species (int)
    return species_nb
end
function getStartSpecie( atoms::T1, specie::T2 ) where { T1 <: AtomList, T2 <: AbstractString }
    start_specie = 0
    nb_atoms = getNbAtoms( atoms )
    for atom=1:nb_atoms
        if specie == atoms.names[atom]
            start_specie = atom
            break
        end
    end
    return start_specie
end
function getStartSpecies( atoms::T1, species::Vector{T2} ) where { T1 <: AtomList, T2 <: AbstractString }
    nb_species = size( species )[1]
    start_species = zeros(Int, nb_species )
    nb_atoms = getNbAtoms( atoms )
    for i_spec=1:nb_species
        for atom=1:nb_atoms
            if species[i_spec] == atoms.names[atom]
                start_species[i_spec] = atom
                break
            end
        end
    end
    return start_species
end
function getStartSpecie( traj::Vector{T1}, specie::T2 ) where { T1 <: AtomList, T2 <: AbstractString }
    nb_step = getNbStep( traj )
    start_specie = zeros(Int, nb_step )
    for step = 1:nb_step
        start_specie[ step ] = getStartSpecie( traj[step], specie )
    end
    return start_specie
end
function getStartSpecies( traj::Vector{T1}, species::Vector{T2} ) where { T1 <: AtomList, T2 <: AbstractString }
    nb_step = getNbStep( traj )
    nb_species = size(species)[1]
    start_species = zeros(Int, nb_step, nb_species  )
    for step=1:nb_step
        start_species[ step, : ] = getStartSpecies( traj[step], species )
    end
    return start_species
end
function getNbElementSpecies( atoms::T1, species::Vector{T2} ) where { T1 <: AtomList, T2 <: AbstractString }
    nb_species = size(species)[1]
    nb_element_species = zeros(Int, nb_species )
    nb_atoms = getNbAtoms( atoms )
    for i_specie = 1:nb_species
        for atom = 1:nb_atoms
            if species[ i_specie ] == atoms.names[atom]
                nb_element_species[ i_specie ] += 1
            end
        end
    end
    return nb_element_species
end
function getNbElementSpecies( traj::Vector{T1}, species::Vector{T2} ) where { T1 <: AtomList, T2 <: AbstractString }
    nb_step = getNbStep( traj )
    nb_species = size(species)[1]
    nb_element_species = zeros( nb_step, nb_species )
    for step = 1:nb_step
        nb_element_species[ step , : ] = getNbElementSpecies( traj[step], species )
    end
    return nb_element_species
end
function getNbElementSpecie( atoms::T1, specie::T2 ) where { T1 <: AtomList, T2 <: AbstractString }
    nb_element_specie = 0
    nb_atoms = getNbAtoms( atoms )
    for atom = 1:nb_atoms
        if  specie == atoms.names[i_specie]
            nb_element_species += 1
        end
    end
    return nb_element_specie
end
function getNbElementSpecie( traj::Vector{T1}, specie::T2 ) where { T1 <: AtomList, T2 <: AbstractString }
    nb_step = getNbStep( traj )
    nb_element_specie = zeros(nb_step)
    for step = 1:nb_step
        nb_element_specie[step]  = getNbElementSpecie( traj[step], specie )
    end
    return nb_element_specie
end
#-------------------------------------------------------------------------------

# Computes velocities by finite different
#-------------------------------------------------------------------------------
function computeVelocities( traj::Vector{T1}, target_step::T2, dt::T3 ) where { T1 <: AtomList, T2 <: Int, T3 <: Real }
    # Argument
    # - traj:
    # - target_step:
    # - dt:

    return ( traj[target_step].positions - traj[target_step-1].positions )/dt
end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
function buildNames( species::Vector{T1}, nb_species::Vector{T2} ) where { T1 <: AbstractString, T2 <: Int }
    nb_atoms=sum(nb_species)
    nb_species_=size(nb_species)[1]
    names=Vector{AbstractString}(undef,nb_atoms)
    for i_spec=1:nb_species_
        start_ = sum( nb_species[1:i_spec-1])
        for atom_spec=1:nb_species[i_spec]
            names[start_+atom_spec] = species[i_spec]
        end
    end
    return names
end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
function makeTrajAtomList( positions::Array{T1,3}, names::Vector{T2}, index::Array{T3} ) where { T1 <: Real, T2 <: AbstractString, T3 <: Int }
    nb_step = size(positions)[1]
    nb_atoms = size(positions)[2]
    traj = Vector{ AtomList }( undef, nb_step)
    for step=1:nb_step
        traj[step] = AtomList( names, index, positions[step,:,:])
    end
    return traj
end
#-------------------------------------------------------------------------------

end
