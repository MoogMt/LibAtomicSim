module topomap

# Load module
using utils

# Description
# - functions used to construct topological maps

# Read topological map frame matrix
#-------------------------------------------------------------------------------
function readFrameToFrameMatrix( file_path::T1 ) where { T1 <: AbstractString }
    # Argument
    # - file_path: path to the input file
    # Output
    # - distance_matrix: matrix of distance between all structures
    # OR false if something is wrong

    # Check that the file exists
    if ! isfile( file_path )
        # If not sends a message and return false, false false
        print("File: ", file_path, " not found!\n")
        return false
    end

    # Opens file
    handle_in = open( file_path )

    # Reads line and parse it with " "
    keywords = split( readline( handle_in ) )

    # Get number of structures
    nb_structure = parse(Int, keywords[1] )

    # Get maximum distance as second element of the first line (parse to Float)
    dist_max = parse(Float64, keywords[2] )

    # Initialize distance matrix
    distance_matrix = zeros(Real, nb_structure, nb_structure )

    # Loop over structures 1
    for structure1=1:nb_structure
        # Prints progress
        print("Reading FRAME_TO_FRAME.MATRIX - progress: ",structure1/nb_structure*100,"%\n" )

        # Read and parse line with " "
        line = split( readline( handle_in ) )

        # Loop over structures 2
        for structure2=1:nb_structure
            # Read and cast distance as float
            distance_matrix[structure1,structure2] = parse(Float64, line[structure2] )*dist_max

            # Use the fact that the matrix is symmetric
            distance_matrix[structure2,structure1] = distance_matrix[structure1,structure2]
        end
    end

    # Close file
    close(handle_in)

    # Returns distance matrix
    return distance_matrix
end
#-------------------------------------------------------------------------------

# Compute cost
#-------------------------------------------------------------------------------
function computeCost( distance_matrix::Array{T1,2}, positions::Array{T2,2}, cost_coeff::T3 ) where { T1 <: Real, T2 <: Real, T3 <: Real }
    # Arguments
    # - distance_matrix: distance_matrix
    # - positions: array with positions in the 2D projection plan
    # - cost_coef: coefficient of the cost the k in k*(x-x0), real
    # Output
    # - cost: errors between actual distances and projection distnaces

    # Get the number of structures
    n_structures = size(distance_matrix)[1]

    # Dimension of the projection
    n_dim = size(positions)[2]

    # Initialize cost at 0
    cost=0.0

    # Loop over structures 1
    for structure1=1:n_structures-1
        # Loop over structure2
        for structure2=structure1+1:n_structures
            # Initialize distance to 0
            dist = 0

            # Loop over dimensions
            for k=1:n_dim
                # Computes distance between positions in the 2D plane for each dimension
                dist_loc = ( positions[i,k] - positions[j,k] )

                # Compute the square and and adds it to the distance
                dist += dist_loc*dist_loc
            end

            # Compute the square root of the distance
            dist = sqrt(dist)

            # Compute distance between the distances between
            # points in projection and points in real space
            dist = dist - distance_matrix[i,j]

            # Add the cost for the couple point
            cost += 0.5*cost_coeff*dist*dist
        end
    end

    # Return the total cost of the projection
    return cost
end
#-------------------------------------------------------------------------------

# Compute errors
#-------------------------------------------------------------------------------
function computeErrors( positions::Array{T1,2}, distances_matrix::Array{T2,2} ) where { T1 <: Real, T2 <: Real }
    # Argument
    # - positions: positions of points in the projection plane
    # - distances_matrix: matrix of actual distances between points
    # Output
    # - mean_err: mean error
    # - max_err: maximum error on the set

    # Get number of points
    nb_points   = size(positions)[1]

    # Get dimension
    n_dimension = size(positions)[2]

    # Initialize max, mean counter
    max_err  = 0
    mean_err = 0

    # Initialize counter
    count_ = 0

    # Loop over structures1
    for structure=1:nb_points-1
        # Loop over structures2
        for structure2=structure+1:nb_points
            # Initialize distance at 0
            distance=0

            # Loop over dimensions
            for i=1:n_dimension
                # Compute distance in dimension i
                dist = positions[ structure, i ] - positions[ structure2, i ]

                # Square distance
                distance += dist*dist
            end

            # Compute square root of distances
            distance = sqrt( distance )

            # Compute error in the projection
            err = abs( distance - distances_matrix[ structure, structure2 ] )

            # Add error to the average
            mean_err += err

            # If err is superior to max value, change the value to max
            if err > max_err
                max_err = err
            end

            # Increments of count
            count_ += 1
        end
    end

    # Returns mean error and max error
    return mean_err/count_, max_err
end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
function monteCarloProjection( n_dim::T1, n_iterations::T2, cost_coeff::T3, move_coef::T4, thermalEnergy::T5, distance_matrix::Array{T6,2} ) where { T1 <: Int,  T2 <: Int, T3 <: Real, T4 <: Real, T5 <: Real, T6 <: Real }
    nb_structure=size(distance_matrix)[1]
    # Randomly put points on the plan
    point_pos=rand(nb_structure,n_dim)
    # Compute initial cost
    cost = computeCost( distance_matrix, point_pos , cost_coeff )
    for iteration=1:n_iterations
        print("MonteCarlo Projection - Progress: ", round(iteration/n_iterations*100,digits=3),"%\n")
        # Choose a random point
        random_point = Int( trunc( rand()*nb_structure+1 ) )
        # Move
        point_pos_moved = copy( point_pos )
        point_pos_moved[ random_point, :] = point_pos[ random_point, : ] .+ (rand(n_dim).-0.5).*move_coef
        # Compute the cost of the move
        cost_move = computeCost( distance_matrix, point_pos_moved , cost_coeff )
        deltaE = ( cost_move - cost )/thermalEnergy
        if deltaE > 0
            # If cost is unfavorable, random
            if rand() < exp(-deltaE)
                point_pos = copy( point_pos_moved )
                cost=cost_move
            end
        else
            # If cost is favorable, accept the move
            point_pos = copy(point_pos_moved )
            cost=cost_move
        end
    end
    return  point_pos, cost
end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
function writeMap( file_out::T1, positions::Array{T2,2} ) where { T1 <: AbstractString, T2 <: Real }
    handle_out = open(file_out,"w")
    nb_structure=size(positions)[1]
    nb_dim = size(positions)[2]
    for structure=1:nb_structure
        for dim=1:nb_dim
            write( handle_out, string( positions[structure,dim], " " ) )
        end
        write( handle_out, string("\n") )
    end
    close( handle_out )
    return true
end
function writeErrors( file_out::T1, positions::Array{T2,2}, distances_matrix::Array{T3,2} ) where { T1 <: AbstractString, T2 <: Real, T3 <: Real }
    handle_out = open( file_out, "w")
    nb_structure = size(positions)[1]
    n_dimension = size(positions)[2]
    for structure=1:nb_structure-1
        for structure2=structure+1:nb_structure
            distance=0
            for i=1:n_dimension
                dist = positions[ structure, i ] - positions[ structure2, i ]
                distance += dist*dist
            end
            distance=sqrt(distance)
            write( handle_out, string( distances_matrix[ structure, structure2 ], " ", distance, "\n" ) )
        end
    end
    close(handle_out)
    return true
end
function writePlotter( file_out::T1, map_file::T2, structure_names::Vector{T3}, positions::Array{T4,2}, columns::Vector{T5}, offset::Vector{T6} ) where { T1 <: AbstractString, T2 <: AbstractString, T3 <: AbstractString, T4 <: Real, T5 <: Int, T6 <: Real }
    nb_point = size(positions)[1]
    handle_out = open( file_out, "w" )
    Base.write( handle_out, string( "set term qt 0 font \"Arial 12,12\" \n" ) )
    for i=1:nb_point
        str = string( "set label \"", structure_names[i] ,"\" at " )
        str = string( str, positions[ i, columns[1] ] + offset[1], "," )
        str = string( str, positions[ i, columns[2] ] + offset[2], "\n" )
        Base.write( handle_out, str)
    end
    Base.write( handle_out, string( "plot \"", map_file, "\" u ", columns[1], ":", columns[2], " ps 1 pt 7 title \"\" \n" ) )
    Base.write( handle_out, string( "set xlabel \"dx\"\n" ) )
    Base.write( handle_out, string( "set ylabel \"dy\"\n" ) )
    close( handle_out )
    return true
end
#-------------------------------------------------------------------------------

end
