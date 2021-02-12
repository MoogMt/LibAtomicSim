module atom_mod

# Loading necessary modules from LibAtomicSim
using utils
using periodicTable

# Export useful functions
# TODO: check that they are all exported
export AtomList, AtomMolList
export switchAtoms!, moveAtom!, sortAtomsByZ!
export getNbAtoms, getNbStep, getNbMol, moveAtom
export getPositionsAsArray

# Description:
# Set of structures and associated function that deals with trajectory in termes of
# atoms, their positions, index and names, and all related functions which help deal
# with all manipulation that does not require or affect the periodic boundary conditions

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
# AtomTraj: Contains a simple trajectory for the atoms
mutable struct AtomTraj

    # Variables
    #-------------------------------
    names::Vector{AbstractString} # Names of the atoms (chemical species)
    index::Vector{Int}            # Index of atoms (labels)
    positions::Array{Real,3}      # Positions of the atoms
    #-------------------------------

    #----------------------------------------------------------------------
    function AtomTraj()
        # Argument
        # - None
        # Output
        # - Creates a trajectory
    end
    #----------------------------------------------------------------------
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

# Create a trajectory as a vector of AtomList from a list of positions, names and index
#-------------------------------------------------------------------------------
function makeTrajAtomList( positions::Array{T1,3}, names::Vector{T2}, index::Array{T3} ) where { T1 <: Real, T2 <: AbstractString, T3 <: Int }
    # Arguments:
    # - positions: Array contains poisitions of atoms (nb_step,nb_atoms,3)
    # - names: vector of string with names of atoms
    # - index: vector of int, with index of atoms
    # Output:
    # traj: vector of AtomList describing the trajectory


    # Get number of steps of trajectory
    nb_step = size(positions)[1]

    # Get number of atoms
    nb_atoms = size(positions)[2]

    # Initialize output
    traj = Vector{ AtomList }( undef, nb_step)

    # Loop over steps
    for step=1:nb_step
        # construct AtomList for each step
        traj[step] = AtomList( names, index, positions[step,:,:])
    end

    # Return trajectory as a vector of AtomList
    return traj
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

    # Return a vector with all species names in the traj, eliminating redundancy
    return unique(species)
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
# Get the first occurence of a target specie
function getStartSpecie( atoms::T1, specie::T2 ) where { T1 <: AtomList, T2 <: AbstractString }
    # Argument
    # - atoms : AtomList containing info on the structure
    # - specie: target specie name
    # Output
    # - start_specie: index of the first occurence of the target specie

    # Initialize output
    start_specie = 0

    # Get the number of atoms
    nb_atoms = getNbAtoms( atoms )

    # Loop over atoms
    for atom=1:nb_atoms
        # If the atom is of the target specie
        if specie == atoms.names[atom]
            # Stock the index of the first occurence and breaks the loop
            start_specie = atom
            break
        end
    end

    # Return the first occurence of the target specie
    return start_specie
end
# Get the first occurence of all species
function getStartSpecies( atoms::T1, species::Vector{T2} ) where { T1 <: AtomList, T2 <: AbstractString }
    # Argument
    # - atoms:   AtomList describing the structure
    # - species: target species
    # Output
    # - start_species: vector of int, with the first occurence of each one of the species

    # Get the number of species
    nb_species = size( species )[1]

    # Initialize output
    start_species = zeros(Int, nb_species )

    # Get the number of atoms
    nb_atoms = getNbAtoms( atoms )

    # Loop over the target species
    for i_spec=1:nb_species
        # Loop over atoms
        for atom=1:nb_atoms
            # if the specie is found
            if species[i_spec] == atoms.names[atom]
                # Stores and break
                start_species[i_spec] = atom
                break
            end
        end
    end

    # Return the first occurences of all species
    return start_species
end
# Get the first occurences of a given specie in a trajectory
function getStartSpecie( traj::Vector{T1}, specie::T2 ) where { T1 <: AtomList, T2 <: AbstractString }
    # Argument
    # - traj: Vector of AtomList, positions
    # - specie: string with name of the target specie
    # Output
    # - start_specie: vector of int, corresponding to the first index occurence of the target specie

    # Get the number of steps
    nb_step = getNbStep( traj )

    # Initialize
    start_specie = zeros(Int, nb_step )

    # Loop over steps
    for step = 1:nb_step
        # For each step stores the first occurence of the target_specie
        start_specie[ step ] = getStartSpecie( traj[step], specie )
    end

    # Return the first occurence of the specie over all the steps of trajectory
    return start_specie
end
# Get the first occurences for each step of several target species for a trajectory
function getStartSpecies( traj::Vector{T1}, species::Vector{T2} ) where { T1 <: AtomList, T2 <: AbstractString }
    # Argument
    # - traj: vector of AtomList as a trajectory
    # - species: target chemical species
    # Output
    # - start_specie: first occurences over the step for all target species

    # Get number of steps of trajectory
    nb_step = getNbStep( traj )

    # Get number of species in the trajectory
    nb_species = size(species)[1]

    # Initialize output
    start_species = zeros(Int, nb_step, nb_species  )

    # Loop over step
    for step=1:nb_step
        # Get first occurence of all target species
        start_species[ step, : ] = getStartSpecies( traj[step], species )
    end

    # Return start_species
    return start_species
end
# Get the number of element of a given specie in an AtomList
function getNbElementSpecie( atoms::T1, specie::T2 ) where { T1 <: AtomList, T2 <: AbstractString }
    # Arguments
    # - atoms: AtomList containing atomic structure
    # - specie: name of the target specie (string)
    # Output
    # -nb_element_specie: Number of elements of the target specie in the structure (int)

    # Initialize the number of specie
    nb_element_specie = 0

    # get the number of atoms
    nb_atoms = getNbAtoms( atoms )

    # Loop over atoms
    for atom = 1:nb_atoms
        # If atom specie is the same as the target...
        if  specie == atoms.names[i_specie]
            # Increment the counter
            nb_element_species += 1
        end
    end

    # Return the counter of number of specie
    return nb_element_specie
end
# Get the number of element for a set of species in an AtomList file
function getNbElementSpecies( atoms::T1, species::Vector{T2} ) where { T1 <: AtomList, T2 <: AbstractString }
    # Argument
    # - atoms: AtomList describing atomic structure
    # - species: target_species
    # Output
    # - nb_element_specie: number of elements for all species

    # Number of different species
    nb_species = size(species)[1]

    # Initialize the number of each species
    nb_element_species = zeros(Int, nb_species )

    # Number of atoms
    nb_atoms = getNbAtoms( atoms )

    # Loop over species
    for i_specie = 1:nb_species
        # Loop over atoms
        for atom = 1:nb_atoms
            # if atom specie corresponds to loop specie ...
            if species[ i_specie ] == atoms.names[atom]
                # Increments the counter of number of elements per specie
                nb_element_species[ i_specie ] += 1
            end
        end
    end

    # returns number of element per species
    return nb_element_species
end
# Get the number of element of a given specie in a traj (vector of AtomList)
function getNbElementSpecie( traj::Vector{T1}, specie::T2 ) where { T1 <: AtomList, T2 <: AbstractString }
    # Argument:
    # - traj: vector of AtomList describing a trajectory
    # - specie: target specie (string)
    # Output
    # nb_element_specie: number of element of the target specie over time in the trajectory

    # Get number of step
    nb_step = getNbStep( traj )

    # Initialize output
    nb_element_specie = zeros(nb_step)

    # Loop over time
    for step = 1:nb_step
        # Get the number of element of the target specie
        nb_element_specie[step]  = getNbElementSpecie( traj[step], specie )
    end

    # Return the number of element of the target specie over time
    return nb_element_specie
end
# Get the number of each element over time in a vector of AtomList
function getNbElementSpecies( traj::Vector{T1}, species::Vector{T2} ) where { T1 <: AtomList, T2 <: AbstractString }
    # Argument
    # - traj: vector of AtomList,
    # - species: vector of string for each target chemical
    # Output
    # - nb_element_species: number of element for each species, over time

    # Get number of steps
    nb_step = getNbStep( traj )

    # Get number of species
    nb_species = size(species)[1]

    # Initialize output
    nb_element_species = zeros( nb_step, nb_species )

    # Loop over steps
    for step = 1:nb_step
        # Get the number of element for each target species for each step
        nb_element_species[ step , : ] = getNbElementSpecies( traj[step], species )
    end

    # return the number of element over time
    return nb_element_species
end
#-------------------------------------------------------------------------------

# Computes velocities by finite different
#-------------------------------------------------------------------------------
function computeVelocities( traj::Vector{T1}, target_step::T2, dt::T3 ) where { T1 <: AtomList, T2 <: Int, T3 <: Real }
    # Argument
    # - traj: Trajectory as a vector of AtomList
    # - target_step: target step for the computation of the velocity
    # - dt: timestep of the trajectory

    # Need at least 2 steps to work
    if target_step < 2
        return false
    end

    # Computes and returns the velocities from finite difference
    return ( traj[target_step].positions - traj[target_step-1].positions )/dt
end
#-------------------------------------------------------------------------------

# Build name vector containing the species names
#-------------------------------------------------------------------------------
function buildNames( species::Vector{T1}, nb_species::Vector{T2} ) where { T1 <: AbstractString, T2 <: Int }
    # Arguments
    # - species: species
    # - nb_species: number of atom for each species
    # Output:
    # - names: names of the species?

    # Get number of atoms in the cell
    nb_atoms=sum(nb_species)

    # Get number of species
    nb_species_=size(nb_species)[1]

    # Initialize output
    names=Vector{AbstractString}(undef,nb_atoms)

    # Loop over species
    for i_spec=1:nb_species_
        # Find where the specie starts in the vector
        start_ = sum( nb_species[1:i_spec-1])
        # Loop over atom
        for atom_spec=1:nb_species[i_spec]
            # Affect all the element of the vector in the range to loop-specie
            names[start_+atom_spec] = species[i_spec]
        end
    end

    # Returns vector of the names of atoms
    return names
end
#-------------------------------------------------------------------------------

# Distance computation (without PBC)
#-------------------------------------------------------------------------------
# Compute distance between two atomic positions (real vectors) without PBC (may be redundant with a function in geom)
function distanceNoPBC( v1::Vector{T1}, v2::Vector{T2} ) where { T1 <: Real, T2 <: Real }
    # Arguments
    # - vector1 : position of atom1
    # - vector2 : position of atom2
    # Output
    # - Distance between the vector 1 and 2

    # Initialize output
    dist=0

    # Loop over dimensions
    for i=1:3
        # Distance in the dimension
        loc = v1[i] - v2[i]
        # add the square of the distance to the total distance
        dist += loc*loc
    end

    # Return the distance between v1 and v2
    return sqrt(dist)
end
# Compute the distance between two atoms (in AtomList frame) without PBC
function distanceNoPBC( atoms::T1, index1::T2, index2::T3 ) where { T1 <: AtomList, T2 <: Int, T3 <: Int }
    # Arguments
    # - atoms: AtomList containing the atomic positions
    # - index1, index2: index of the target atoms
    # Output
    # - distance between atoms index1 and index2

    # Computes and return the distance betweens the atoms
    return distanceNoPBC( atoms.positions[index1], atoms.positions[index2,:] )
end
#-------------------------------------------------------------------------------

end
