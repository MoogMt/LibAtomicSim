module topomap

using utils

#-------------------------------------------------------------------------------
function computeCost( distance_matrix::Array{T1,2}, positions::Array{T2}, cost_coeff::T3 ) where { T1 <: Real, T2 <: Real, T3 <: Real }
    n_structures=size(distance_matrix)[1]
    cost=0.
    n_dim=size(positions)[2]
    for i=1:n_structures-1
        for j=i+1:n_structures
            dist = 0
            for k=1:n_dim
                dist += ( positions[i,k] - positions[j,k] )*( positions[i,k] - positions[j,k] )
            end
            dist=sqrt(dist)
            cost += 0.5*cost_coeff*( dist - distance_matrix[i,j] )*( dist - distance_matrix[i,j] )
        end
    end
    return cost
end
function readFrameToFrameMatrix( file_path::T1 ) where { T1 <: AbstractString }
    if ! isfile(file_path)
        print("File: ", file_path, " not found!\n")
        return false, false, false
    end
    handle_in = open( file_path )
    keywords=split( readline( handle_in ) )
    nb_structure = parse(Int, keywords[1] )
    dist_max = parse(Float64, keywords[2] )
    distance_matrix=zeros(nb_structure,nb_structure)
    for i=1:nb_structure
        print("Reading FRAME_TO_FRAME.MATRIX - progress: ",i/nb_structure*100,"%\n" )
        line= split( readline( handle_in ) )
        for j=1:nb_structure
            distance_matrix[i,j] = parse(Float64, line[j] )
            distance_matrix[j,i] = distance_matrix[i,j]
        end
    end
    close(handle_in)
    return distance_matrix, dist_max, nb_structure
end
function monteCarloProject( n_dim::T1, n_iterations::T2, cost_coeff::T3, move_coef::T4, thermalEnergy::T5, distance_matrix::Array{T6,2} ) where { T1 <: Int,  T2 <: Int, T3 <: Real, T4 <: Real, T5 <: Real, T6 <: Real }
    nb_structure=size(distance_matrix)[1]
    # Randomly put points on the plan
    point_pos=rand(nb_structure,n_dim)
    # Compute initial cost
    cost = computeCost( distance_matrix, point_pos , cost_coeff )
    for iteration=1:n_iterations
        print("MonteCarlo Projection - Progress: ",iteration/n_iterations*100,"%\n")
        # Choose a random point
        random_point = round(Int, rand()*nb_structure+1 ) #
        # Move
        point_pos_moved = copy( point_pos )
        point_pos_moved[ random_point, :] = point_pos[ random_point, : ] .+ (rand(n_dim).-0.5).*move_coef
        # Compute the cost of the move
        cost_move = computeCost( distance_matrix, point_pos_moved , cost_coeff )
        deltaE = ( cost_move - cost )/thermalEnergy
        if de > 0
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
    for structure=1:nb_structure
        for structure2=structure+1:nb_structure
            distance=0
            for i=1:n_dimension
                distance += (positions[structure,i]-positions[structure2,i])*(positions[structure,i]-positions[structure2,i])
            end
            distance=sqrt(distance)
            write( file_out, string( distances_matrix[structure,structure2], " ", distance, "\n" ) )
        end
    end
    close(handle_out)
    return true
end
#-------------------------------------------------------------------------------

end
