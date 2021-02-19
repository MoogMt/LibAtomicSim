module topomap

# Load module
using utils

# Description
# - functions used to construct topological maps

# Exporting functions
export readFrameToFrameMatrix
export computeCost
export computeErrors
export monteCarloProjection
export writeMap, writeErrors, writePlotter

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

# Use Monte Carlo to project points in n_dimension space
#-------------------------------------------------------------------------------
function monteCarloProjection( n_dim::T1, n_iterations::T2, cost_coeff::T3, move_coef::T4, thermalEnergy::T5, distance_matrix::Array{T6,2} ) where { T1 <: Int,  T2 <: Int, T3 <: Real, T4 <: Real, T5 <: Real, T6 <: Real }
    # Arguments
    # - n_dim: number of dimension (int)
    # - n_iterations: maximum number of iterations of Monte Carlo
    # - cost_coef: coefficient for the cost of error
    # - move_coef : coefficient that decides the largeness of the Monte Carlo move
    # - thermalEnergy: thermal energy of the monte carlo (real,scalar)
    # - distance_matrix: matrix of the actual distances (real,nb_points,nb_points)
    # Output:
    # - point_pos: posiitons of points in the projections
    # - cost: final cost of the projection

    # Get number of structures
    nb_structure = size(distance_matrix)[1]

    # Randomly put points on the plan
    point_pos=rand(nb_structure,n_dim)

    # Compute initial cost
    cost = computeCost( distance_matrix, point_pos , cost_coeff )

    # Loop over iterations
    for iteration=1:n_iterations
        # Print progression of the Monte-Carlo projection
        print("MonteCarlo Projection - Progress: ", round(iteration/n_iterations*100,digits=3),"%\n")

        # Choose a random point
        random_point = Int( trunc( rand()*nb_structure+1 ) )

        # Move
        point_pos_moved = copy( point_pos )
        point_pos_moved[ random_point, :] = point_pos[ random_point, : ] .+ (rand(n_dim).-0.5).*move_coef

        # Compute the cost of the move
        cost_move = computeCost( distance_matrix, point_pos_moved , cost_coeff )

        # Compute difference in energy
        deltaE = ( cost_move - cost )/thermalEnergy

        # If the deltaE is superior to 0 then Monte Carlo choice
        if deltaE > 0
            # If cost is unfavorable, randomly chose to accept it or not
            if rand() < exp(-deltaE)
                # If we accept it, we copy the position and cost
                point_pos = copy( point_pos_moved )
                cost = cost_move
            end
        # If cost is favorable, accept the move
        else
            # Copy the postiions and the cost
            point_pos = copy(point_pos_moved )
            cost = cost_move
        end
    end

    # positions of the points in the projection, cost of the projection
    return  point_pos, cost
end
#-------------------------------------------------------------------------------

# Writing data to file
#-------------------------------------------------------------------------------
# Writting positions of the projection to file
function writeMap( file_path::T1, positions::Array{T2,2} ) where { T1 <: AbstractString, T2 <: Real }
    # Argument
    # - file_path: path to the output file
    # - positions: positions of the points in the projection plan
    # Output
    # - Bool: whether the thing was right

    # Opens output file
    handle_out = open( file_out, "w" )

    # Get the nb of structures
    nb_structure=size(positions)[1]

    # Get number of dimension of the projection
    nb_dim = size(positions)[2]

    # Loop over the structures
    for structure=1:nb_structure
        # Loop over dimensions
        for dim=1:nb_dim
            # Write positions to file
            write( handle_out, string( positions[structure,dim], " " ) )
        end
        # Write end of line to file
        write( handle_out, string("\n") )
    end

    # Closing output file
    close( handle_out )

    # Return true if everything is right
    return true
end
# Writting errors to file
function writeErrors( file_path::T1, positions::Array{T2,2}, distances_matrix::Array{T3,2} ) where { T1 <: AbstractString, T2 <: Real, T3 <: Real }
    # Argument
    # - file_path: path of the output file
    # - positions: positions of the points in projections
    # - distances_matrix: matrix of the distances between points in actual space
    #  Output
    # - Bool, whether the writting went ok

    # Get number of structures
    nb_structure = size(positions)[1]

    # Get number of dimensions
    n_dimension = size(positions)[2]

    # Open output file
    handle_out = open( file_path, "w" )

    # Loop over structure1
    for structure=1:nb_structure-1
        # Loop over structure2
        for structure2=structure+1:nb_structure
            # Initialize distances to 0
            distance=0

            # Loop over dimensions
            for i=1:n_dimension
                # Compute distances in i-th dimension
                dist = positions[ structure, i ] - positions[ structure2, i ]

                # Adds square of the distances in i-th dimension
                distance += dist*dist
            end

            # Compute square root of the distance
            distance = sqrt( distance )

            # Write distances in actual and projection space
            write( handle_out, string( distances_matrix[ structure, structure2 ], " ", distance, "\n" ) )
        end
    end

    # Close ouput file
    close( handle_out )

    # Return true if we managed to get to that point
    return true
end
# Writting plotter file to plot map using Gnuplot
function writePlotter( file_out::T1, map_file::T2, point_names::Vector{T3}, positions::Array{T4,2}, columns::Vector{T5}, offset::Vector{T6} ) where { T1 <: AbstractString, T2 <: AbstractString, T3 <: AbstractString, T4 <: Real, T5 <: Int, T6 <: Real }
    # Argument
    # - file_path: path of the output file
    # - map_file: path to the file that we want to load in gnuplot
    # - point_names: name of each point
    # - positions: positions in the projection space of the points
    # - columns: index of the colums that will be used to plot data
    # - offset: offset in (x,y) dimension of the text
    # Output
    # - Bool, wether the writing was ok

    # Get number of points
    nb_point = size( positions )[1]

    # Opens output file
    handle_out = open( file_out, "w" )

    # Write setting line
    Base.write( handle_out, string( "set term qt 0 font \"Arial 12,12\" \n" ) )

    # Loop over points
    for point=1:nb_point
        # Create string for priting each label in space near to the points
        str = string( "set label \"", structure_names[point] ,"\" at " )
        str = string( str, positions[ point, columns[1] ] + offset[1], "," )
        str = string( str, positions[ point, columns[2] ] + offset[2], "\n" )

        # Write string to file
        Base.write( handle_out, str)
    end

    # Write plotting line
    Base.write( handle_out, string( "plot \"", map_file, "\" u ", columns[1], ":", columns[2], " ps 1 pt 7 title \"\" \n" ) )

    # Write labels on x and y labels
    Base.write( handle_out, string( "set xlabel \"dx\"\n" ) )
    Base.write( handle_out, string( "set ylabel \"dy\"\n" ) )

    # Closing output file
    close( handle_out )

    # Returns true if all went well
    return true
end
#-------------------------------------------------------------------------------

end
