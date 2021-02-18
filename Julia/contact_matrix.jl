module contact_matrix

# Loading necessary modules
using atom_mod
using cell_mod
using graph

# Description
# This functions allows the calculation of distance matrix (containing all distances
# between all atoms in the structure). But also the computation of adjacency/contact
# matrix that contains information about the bonding/proximity between all atoms.
#
# The different between the three types being primarly that a logistic function is used
# to get the proximity and normalize the distances of the distace matrix in a contact matrix
# while a heavyside function is used to get the adjacency matrix

# Export functions

export computeDistanceMatrix, computeAdjacencyMatrix
export readMatrix
export writeMatrix

# Building Distance Matrix
# - contains distances between all atoms/points in the structure/set
#-------------------------------------------------------------------------------
# Compute the distance matrix for an atomic structure
function computeDistanceMatrix( atoms::T1, cell::T2 ) where { T1 <: atom_mod.AtomList , T2 <: cell_mod.Cell_param }
    # Arguments
    # - atoms: AtomList containing the atomic positions and (names and index)
    # - cell: Cell_param containing cell informations
    # Output
    # - distance_matrix: real array (nb_atoms,nb_atoms), contains the distance
    # matrix

    # Get the number of atoms in the structure
    nb_atoms = size(atoms.names)[1]

    # Initialize distance matrix
    distance_matrix = zeros(Real, nb_atoms, nb_atoms )

    # Loop over atoms 1
    for atom1=1:nb_atoms
        # Loop over atoms 2
        for atom2=atom1+1:nb_atoms
            # Compute distances between atom1 and atom2
            distance_matrix[atom1,atom2] = cell_mod.distance( atoms, cell, atom1, atom2 )

            # The matrix is symmetric so we have M[i,j] = M[j,i]
            distance_matrix[atom2,atom1] = distance_matrix[atom1,atom2]
        end
    end

    # Returns the distance matrix
    return distance_matrix
end
# Compute the distance matrix for a set of points in arbitrary space
function computeDistanceMatrix( data::Array{T1} ) where { T1 <: Real }
	# Argument
	# - data: Array (real, (nb_point,n_dim) )
	# Output
	# - distance_matrix: real matrix (nb_point,nb_point) containing distances
	# between all points

	# Get the number of data points in the set
	size_data = size(data)[1]

	# Initialize distance matrix
    distance_matrix = zeros(Real, size_data, size_data )

	# Loop over data points 1
	for point1=1:size_data
		# Loop over data point 2
        for point2=point1+1:size_data
			# Computes the distances between point1 and point2
            distance_matrix[point1,point2] = computeDistance( data, point1, point2 )

			# Use the fac that the matrix is symmetric
			distance_matrix[point1,point2] = distance_matrix[point2,point1]
        end
    end

	# Returns distance_matrix
    return distance_matrix
end
# Compute the distance matrix for a set of points in vector space normalizing the distances using min and max of data
function computeDistanceMatrix( data::Array{T1}, max::T2, min::T3 ) where { T1 <: Real, T2 <: Real, T3 <: Real }
	# Argument
	# - data: Vector containing the data
	# - max, min: maximum and minumum values of the data
	# Output
	# - distance_matrix: matrix containing the normalized distances

	# get number of data points
	size_data = size(data)[1]

	# Initialize distance matrix
    distance_matrix = zeros(Real, size_data, size_data )

	# Compute delta between max and min
	delta = max - min

	# Copy and normalize data
	data_copy = copy(data)
	data_copy = (data .- min )./delta

	# Loop over data points 1
    for point1=1:size_data
		# Loop over data points 2
        for point2=1:size_data
			# Compute distance between point1 and point2 normalized
			distance_matrix[point1,point2] = norm( data_copy[point1] - data_copy[point2] )

			# Use the fact that the matrix is symmetric
			distance_matrix[point2,point1] = distance_matrix[point1,point2]
        end
    end

	# Retun distance_matrix
    return distance_matrix
end
# Compute the distance matrix for a set of points in arbitrary space, normalizing the distances using min and max of data
function computeDistanceMatrix( data::Array{T1,2}, max::Vector{T3}, min::Vector{T4} ) where { T1 <: Real, T2 <: Int, T3 <: Real, T4 <: Real }
	# Argument
	# - data:
	# - max, min: maximum and minimum value for normalization
	# Output
	# - matrix: real matrix containg normalized distances

	# Get the size of the data set
	size_data = size(data)[1]

	# Initialize distance matrix
    distance_matrix = zeros(Real, size_data, size_data )

	# Get delta between max and min
	delta = max-min

	# Copy and normalize data
	data_copy = copy(data)
	data_copy = ( data .- min )./delta

	# Loop over data points 1
    for i=1:size_data
		# Loop over data points 2
        for j=1:size_data
			# Compute distances between point 1 and 2
            distance_matrix[i,j] = norm( data_copy[point1,:] .- data_copy[point2,:] )

			# Use the fact that the matrix is symmetric
			distance_matrix[j,i] = distance_matrix[i,j]
        end
    end

	# Returns distance matrix
    return distance_matrix
end
#-------------------------------------------------------------------------------

# Computes Adjacency matrix
# - Contains information whether or not atoms are bonded or not
#-------------------------------------------------------------------------------
function computeAdjacencyMatrix( atoms::T1 , cell::T2, cut_off::T3 ) where { T1 <: atom_mod.AtomList , T2 <: cell_mod.Cell_param, T3 <: Real }
    # Argument
    # - atoms: AtomList with atomic positions
    # - cell: Cell_param with cell informations
    # - cut_off: (real scalar) cut-off to determine bonding between atoms
    # Output
    # - adjacency_matrix: matrix containing whether or not atoms are bonded

    # Get number of atoms in the structure
    nb_atoms = size(atoms.names)[1]

    # Initialize Adjacency matrix
    adjacency_matrix = zeros(Int, nb_atoms, nb_atoms )

    # Loop over atoms 1
    for atom1=1:nb_atoms
        # Loop over atoms 2
        for atom2=i+1:nb_atoms
            # Computes the distance between atom 1 and 2
            # If it is lower than the cut-off
            # Atoms are bonded, and matrix is put to 1
            if cell_mod.distance( atoms, cell, i, j ) <= cut_off
                adjacency_matrix[i,j] = 1
                adjacency_matrix[j,i] = 1
            end
        end
    end

    # Returns adjacency matrix
    return adjacency_matrix
end
#-------------------------------------------------------------------------------

# Read Distance matrix files
#-------------------------------------------------------------------------------
function readStepMatrix( handle_in::T1 , nb_atoms::T2) where { T1 <: IO, T2 <: Int }
    matrix=zeros(nb_atoms,nb_atoms)
    for i=1:nb_atoms
        keywords = split( readline( handle_in ) )
        for j=1:nb_atoms
            matrix[i,j] = parse( Float64, keywords[j] )
        end
    end
    return matrix
end
function readMatrix( file::T1, target_step::T2 ) where { T1 <: AbstractString, T2 <: Int }
    handle_in = open( file )
    # Getting info
    keywords = split( readline( file ) )
    nb_step  = parse(Int, keywords[1] )
    nb_atoms = parse(Int, keywords[2] )
    if step > nb_step
        print("The trajectory is only ", nb_step, " long, you're asking for step ", target_step,".\n")
        print("Stopping now!\n")
        return false
    end
    # Reading
    matrix=zeros(nb_atoms,nb_atoms)
    for step=1:nb_step
        if step == target_step
            matrix[ :, : ]  = readStepMatrix( handle_in, nb_atoms )
        else
            readStepMatrix( handle_in, nb_atoms )
        end
    end
    close(file)
    return matrix
end
function readMatrix( file::T1 ) where { T1<: AbstractString }
    handle_in = open( file )
    # Getting info
    keywords = split( readline( file ) )
    nb_step  = parse(Int, keywords[1] )
    nb_atoms = parse(Int, keywords[2] )
    # Reading
    matrix=zeros(nb_step,nb_atoms,nb_atoms)
    for step=1:nb_step
        matrix[ step, :, : ]  = readStepMatrix( handle_in, nb_atoms )
    end
    close(file)
    return matrix
end
#-------------------------------------------------------------------------------

# Writting distance matrix files
#-------------------------------------------------------------------------------
function writeStepMatrix( handle_out::T1, matrix::Array{T2,2} ) where { T1 <: IO , T2 <: Real }
    nb_atoms=size(matrix)[1]
    for atom1=1:nb_atoms
        for atom2=1:nb_atoms
            write( handle_out, string( matrix[ atom1, atom2 ], " " ) )
        end
        write( handle_out , "\n")
    end
    return true
end
function writeMatrix( file_out::T1, matrix::Array{T2,3} ) where { T1 <: AbstractString, T2 <: Real }
    handle_out = open( file_out, "w" )
    nb_step = size( matrix )[1]
    nb_atoms = size( matrix )[2]
    Base.write( handle_out, string( nb_step, " ", nb_atoms, "\n" ) )
    for step=1:nb_atoms
        writeStepMatrix( handle_out, matrix[step,:,:] )
    end
    close( handle_out )
    return true
end
#-------------------------------------------------------------------------------

end
