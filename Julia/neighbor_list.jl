module neighbor_list

# Loading necessary modules from LibAtomicSim
using utils
using geom
using periodicTable

# Description:
# Set of structures and associated function that deals with trajectory in termes of
# atoms, their positions, index and names, and all related functions which help deal
# with all manipulation that does not require or affect the periodic boundary conditions

# Structures
#-------------------------------------------------------------------------------
# A neighbour list
mutable struct NeighborList

    # Variables
    #-------------------------------
    lut_neighbors_i::Vector{Int} # Names of the atoms (chemical species)
    lut_neighbors_j::Vector{Int}
    number_neighbors::Vector{Int}
    cum_number_neighbors::Vector{Int}
    #-------------------------------

    # Constructors
    #----------------------------------------------------------------------------
    # Creates a default NeighborList
    function NeighborList()
        # Arguments
        # None
        # Output
        # - Creates a default AtomList

        # Create the new AtomList with nothing in it
        new( zeros(0),   # LUT for first neighbor
             zeros(0),   # LUT for first neighbor
             zeros(0),   # Number of neighbors for each point
             zeros(0) )  # Cumulative number of neighbors
    end
    # Creates a default NeighborList
    function NeighborList()
        # Arguments
        # None
        # Output
        # - Creates a default AtomList

        # Create the new AtomList with nothing in it
        new( zeros(0),   # LUT for first neighbor
             zeros(0),   # LUT for first neighbor
             zeros(0),   # Number of neighbors for each point
             zeros(0) )  # Cumulative number of neighbors
    end
    #----------------------------------------------------------------------------
end

end
