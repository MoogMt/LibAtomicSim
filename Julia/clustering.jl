module clustering

export computeDistance, computeDistanceMatrix, computeDistanceMatrixAndMax, swap
export initializeCenter, voronoiAssign, voronoiAssignSingle, voronoiAssignAll
export computeCost, updateCenters, sortMatrix
export kmedoidClustering, computeClusteringCoefficients
export dauraClustering
export gaussianKernel, maxArray, simpleSequence, maxVector
export densityPeakClusteringTrain
export createBlob

function createBlob( n_point::T1 ) where { T1 <: Int }
	points=zeros(n_point)
	return points
end

function createBlobs( n_points::Vector{T1} , centers::Array{T2,2}, spread::Vector{T3} ) where { T1 <: Int, T2 <: Real, T3 <: Real }
	n_blobs=size(n_points)[1]
	n_points_total=sum(n_points)
	n_dim = size(centers)[2]
	points=zeros(sum(n_points_total),n_dim)
	spread_2 = spread.*spread
	for i=1:n_blobs
		start_count=sum(n_points[1:i-1])
		for j=1:n_points[i]
			try_ = (rand(n_dim).-0.5)*2*spread[i]
			while sum(try_.*try_) > spread_2[i]
				try_ = (rand(n_dim).-0.5)*2*spread[i]
			end
			points[start_count+j,:] = centers[i,:] .+ try_
		end
	end
	return points
end

function createRing( n_points::T1, centers::Vector{T2}, small_radius::T3, width::T4 ) where { T1 <: Int, T2 <: Real, T3 <: Real, T4 <: Real }
	n_dim=size(centers)[1]
	points=zeros(Real, n_points,n_dim)
	R=0
	angles=ones(Real, n_dim-1)
	for i=1:n_points
		# Randomize a distnace to the center
		R = small_radius+rand()*width
		angles=rand(n_dim-1)*pi
		angles[n_dim-1]=rand()*2*pi
		for j=1:n_dim-1
			points[i,j]=R*cos(angles[j])
		end
		points[i,n_dim]=R
		for j=2:n_dim
			for k=1:j-1
				points[i,j] = points[i,j]*sin(angles[k])
			end
		end
	end
	return points
end

function computeDistance( data::Array{T1}, data_point::Vector{T2}, index::T3 ) where { T1 <: Real, T2 <: Real, T3 <: Int }
	return sum( ( data[ index, :  ] - data_point[:] ).*( data[ index, :  ] - data_point[:] ) )
end
function computeDistance( data::Array{T1}, index1::T2, index2::T3 ) where { T1 <: Real, T2 <: Int, T3 <: Int }
	return sum( ( data[ index1, : ]-data[ index2 , : ] ).*( data[ index1, : ]-data[ index2 , : ] ) )
end
function computeDistance( data::Array{T1}, n_dim::T2 , max::Vector{T3}, min::Vector{T4}, index1::T5, index2::T6 ) where { T1 <: Real,T2 <: Int,  T3 <: Real, T4 <: Real, T5 <: Int, T6 <: Int }
    dist=0
    for i=1:n_dim
        dist += ( ((data[index1,i]-min[i])/(max[i]-min[i])) - ((data[index2,i]-min[i])/(max[i]-min[i])) )*( ((data[index1,i]-min[i])/(max[i]-min[i])) - ((data[index2,i]-min[i])/(max[i]-min[i])) )
    end
	return dist
end
function computeDistanceMatrix( data::Array{T1} ) where { T1 <: Real }
    size_data=size(data)[1]
    matrix=zeros(size_data,size_data)
    for i=1:size_data
        for j=1:size_data
            matrix[i,j] = computeDistance( data, i, j )
        end
    end
    return matrix
end
function computeDistanceMatrix( data::Array{T1}, max::T2, min::T3 ) where { T1 <: Real, T2 <: Real, T3 <: Real }
    size_data=size(data)[1]
    matrix=zeros(size_data,size_data)
    for i=1:size_data
        for j=1:size_data
			matrix[i,j] = sqrt( ( ((data[i]-min)/(max-min)) - ((data[j]-min)/(max-min)) )*( ((data[i]-min)/(max-min)) - ((data[j]-min)/(max-min)) ) )
        end
    end
    return matrix
end
function computeDistanceMatrix( data::Array{T1}, n_dim::T2, max::Vector{T3}, min::Vector{T4} ) where { T1 <: Real, T2 <: Int, T3 <: Real, T4 <: Real }
    size_data=size(data)[1]
    matrix=zeros(size_data,size_data)
    for i=1:size_data
        for j=1:size_data
            matrix[i,j] = computeDistance( data, n_dim, max, min, i, j )
        end
    end
    return matrix
end
function computeDistanceMatrixAndMax( data::Array{T1}, n_dim::T2, max::Vector{T3}, min::Vector{T4} ) where { T1 <: Real, T2 <: Int, T3 <: Real, T4 <: Real }
    size_data=size(data)[1]
	max_matrix=0
    matrix=zeros(size_data,size_data)
    for i=1:size_data
        for j=1:size_data
            matrix[i,j] = computeDistance( data, n_dim, max, min, i, j )
			if matrix[i,j] > max_matrix
				max_matrix = matrix[i,j]
			end
        end
    end
    return matrix, max_matrix
end
function swap( table::Vector{T1}, index1::T2, index2::T3 ) where { T1 <: Real, T2 <: Int, T3 <: Int }
    stock=table[index1]
    table[index1]=table[index2]
    table[index2]=stock
    return
end
function initializeCenters( n_structures::T1 , distance_matrix::Array{T2,2}  , n_clusters::T3  ) where { T1 <: Int, T2 <: Real , T3 <: Int}
    # Bookkeep
    available=ones(Int,n_structures)

    cluster_centers=zeros(Int,n_clusters)
    prob=zeros(n_structures)

    cluster_centers[1]= trunc( rand()*n_structures )   + 1
    available[ cluster_centers[1] ] = 0
    for i=2:n_clusters
        for k=1:n_structures
            min_distance=1
            for j=1:i-1
                dist_temp=distance_matrix[ k, Int(cluster_centers[j]) ]
                if dist_temp < min_distance
                    min_distance = dist_temp
                end
            end
            prob[k] = min_distance^2
        end
        prob=cumsum(prob/sum(prob))
        found=false
        while ! found
            r=rand()
            if r < prob[1] && available[1] == 1
                cluster_centers[i] = 1
                available[1] = 0
                found = true
            else
                for k=2:n_structures
                    if r > prob[k-1] && r < prob[k] && available[k] == 1
                        cluster_centers[i] = k
                        available[k] = 0
                        found = true
                        break
                    end
                end
            end
        end
    end

    return cluster_centers
end
# Voronoi assignment of points
function voronoiAssign( data::Array{T1}, n_clusters::T2 , cluster_centers::Vector{T3}, data_points::Array{T4} ) where { T1 <: Real, T2 <: Int, T3 <: Int, T4 <: Real , T5 <: Real, T6 <: Real }
	nb_data=size(data)[1]
	dim_data=size(data)[2]
	max=zeros(dim_data)
	for i=1:dim_data
		for j=1:nb_data
			if max[i] < data[j,i]
				max[i] = data[j,i]
			end
		end
	end
	min=max
	for i=1:dim_data
		for j=1:nb_data
			if min[i] > data[j,i]
				min[i] = data[j,i]
			end
		end
	end
	n_points=size(data_points)[1]
	point_clusters=zeros(Int,n_points)
	for i=1:n_points
		point_clusters[i] = voronoiAssign( data, n_clusters, cluster_centers, data_points[i,:], max, min )
	end
	return point_clusters
end
function voronoiAssign( data::Array{T1}, n_clusters::T2 , cluster_centers::Vector{T3}, data_point::Vector{T4} , max::Vector{T5}, min::Vector{T6} ) where { T1 <: Real, T2 <: Int, T3 <: Int, T4 <: Real , T5 <: Real, T6 <: Real }
	index_cluster=1
	min_dist=sum( ( (data[ cluster_centers[1], :  ]-min)./(max-min) - (data_point-min)./(max-min) ).*( (data[ cluster_centers[1], :  ]-min)./(max-min) - (data_point-min)./(max-min) ) )
	for i=2:n_clusters
		dist=sum( ( (data[ cluster_centers[i], :  ]-min)./(max-min) - (data_point-min)./(max-min) ).*( (data[ cluster_centers[i], :  ]-min)./(max-min) - (data_point-min)./(max-min) ) )
		if dist < min_dist
			min_dist=dist
			index_cluster=i
		end
	end
	return index_cluster
end
function voronoiAssign( data::Array{T1}, n_clusters::T2 , cluster_centers::Vector{T3}, data_points::Array{T4} , max::Vector{T5}, min::Vector{T6} ) where { T1 <: Real, T2 <: Int, T3 <: Int, T4 <: Real , T5 <: Real, T6 <: Real }
	n_points=size(data_points)[1]
	point_clusters=zeros(Int,n_points)
	for i=1:n_points
		point_clusters[i] = voronoiAssign( data, n_clusters, cluster_centers, data_points[i,:], max, min )
	end
	return point_clusters
end
function voronoiAssignSingle( distance_matrix::Array{T1,2}, nb_clusters::T2, cluster_centers::Vector{T3}, index::T4 ) where { T1 <: Real, T2 <: Int, T3 <: Int, T4 <: Int}
    min_dist=distance_matrix[ index ,cluster_centers[1] ]
    index_cluster=1
    for i=2:nb_clusters
    dist=distance_matrix[ index ,cluster_centers[i] ]
        if dist < min_dist
            min_dist = dist
            index_cluster = i
        end
    end
    return index_cluster
end
function voronoiAssignAll( nb_structures::T1 , distance_matrix::Array{T2,2}, nb_clusters::T3, cluster_centers::Vector{T4} ) where { T1 <: Int, T2 <: Real, T3 <: Int, T4 <: Int}
    # Index of the cluster for each structure
    cluster_indexs = zeros(Int, nb_structures )
    assignments = zeros(Int, nb_clusters, nb_structures )
    # Contains sizer of clusters
    cluster_sizes=zeros(Int,nb_clusters)
    for structure=1:nb_structures
        cluster_indexs[ structure ] = voronoiAssignSingle( distance_matrix, nb_clusters, cluster_centers, structure )
        cluster_sizes[ cluster_indexs[structure] ] += 1
        assignments[ cluster_indexs[structure], cluster_sizes[ cluster_indexs[structure] ] ] = structure
    end
    return cluster_indexs, cluster_sizes, assignments
end
function voronoiAssignAll( n_structures::T1 , distance_matrix::Array{T2,2}, nb_clusters::T3, cluster_centers::Vector{T4}, cluster_indexs::Vector{T5}, cluster_sizes::Vector{T6}, assignments::Array{T7,2} ) where { T1 <: Int, T2 <: Real, T3 <: Int, T4 <: Int, T5 <: Int, T6 <: Int, T7 <: Int }
    cluster_sizes=zeros(Int,nb_clusters)
    for structure=1:n_structures
        cluster_indexs[ structure ] = voronoiAssignSingle( distance_matrix, nb_clusters, cluster_centers, structure )
        cluster_sizes[ cluster_indexs[structure] ] += 1
        assignments[ cluster_indexs[ structure ], cluster_sizes[ cluster_indexs[ structure ] ] ] = structure
    end
    return
end
function computeCost( n_structures::T1, distance_matrix::Array{T2,2}, cluster_centers::Vector{T3} , cluster_indexs::Vector{T4} ) where { T1 <: Int, T2 <: Real, T3 <: Int, T4 <: Int }
    cost=0
    for i=1:n_structures
        cost += distance_matrix[ i, cluster_centers[ cluster_indexs[i] ] ]
    end
    return cost
end
function updateCenters( distance_matrix::Array{T1,2} , n_clusters::T2, cluster_centers::Vector{T3}, cluster_sizes::Vector{T4}, assignments::Array{T5,2} ) where { T1 <: Real, T2 <: Int, T3 <: Int, T4 <: Int , T5 <: Int }
    for cluster=1:n_clusters
        new_center = cluster_centers[ cluster ]
        cost_min = sum(distance_matrix[ assignments[ cluster,1], assignments[ cluster ,1:cluster_sizes[ cluster ] ] ])
        for structure_in_cluster=2:cluster_sizes[ cluster ]
            cost = sum(distance_matrix[ assignments[ cluster, structure_in_cluster ], assignments[ cluster ,1:cluster_sizes[ cluster ] ] ])
            if cost < cost_min
                cost_min=cost
                new_center = assignments[ cluster , structure_in_cluster ]
            end
        end
        cluster_centers[ cluster ] = new_center
    end
end


# K-menoid
#============================================================================#
# algorithm from H.S.Park and C.H.Jun, Expert Syst. Appl. 36, 3336, 2009.
function kmedoidClustering( n_structures::T1 , distance_matrix::Array{T2,2}, n_clusters::T3 ) where { T1 <: Int, T2 <: Real, T3 <: Int }
    # Initialization of centers
    cluster_centers=initializeCenters(n_structures, distance_matrix, n_clusters )
    # Assign all clusters
    cluster_indexs, cluster_sizes, assignments =voronoiAssignAll( n_structures, distance_matrix, n_clusters, cluster_centers )
    # Compute original cost
    old_cost=computeCost(n_structures,distance_matrix,cluster_centers,cluster_indexs)

    old_cost=1
    while true
        voronoiAssignAll( n_structures, distance_matrix, n_clusters, cluster_centers,cluster_indexs, cluster_sizes, assignments  )
        cost=computeCost( n_structures, distance_matrix, cluster_centers, cluster_indexs)
		if cost < old_cost
            break
        end
        old_cost=cost
        updateCenters( distance_matrix, n_clusters, cluster_centers, cluster_sizes, assignments )
    end
    return cluster_indexs, cluster_centers, cluster_sizes
end
function kmedoidClustering( n_structures::T1 , distance_matrix::Array{T2,2}, n_clusters::T3, n_repeat::T4 ) where { T1 <: Int, T2 <: Real, T3 <: Int, T4 <: Int }
    # First Kmenoid
    cluster_indexs_best, cluster_centers_best, cluster_sizes_best = kmedoidClustering( n_structures, distance_matrix, n_clusters )
    old_cost=computeCost( n_structures, distance_matrix, cluster_centers_best, cluster_indexs_best )
    for i=2:n_repeat
        cluster_indexs, cluster_centers, cluster_sizes = kmedoidClustering( n_structures, distance_matrix, n_clusters )
        cost=computeCost( n_structures, distance_matrix, cluster_centers, cluster_indexs )
        if old_cost > cost
            cluster_indexs_best = cluster_indexs
            cluster_centers_best = cluster_centers
            cluster_sizes_best = cluster_sizes
        end
    end
    return cluster_indexs_best, cluster_centers_best, cluster_sizes_best
end
#==============================================================================#
# Inspired by PIV_clustering
# by G.A. Gallet and F. Pietrucci, 2014
# J.Chem.Phys., 139 , 074101, 2013
function computeClusteringCoefficients( distance_matrix::Array{T1,2}, n_clusters::T2 , cluster_sizes::Vector{T3} , assignments::Array{T4,2} , cut_off::T5 ) where { T1 <: Real, T2 <: Int, T3 <: Int, T4 <: Int , T5 <: Real }
    clustering_coefficients=zeros(n_clusters)
    if cut_off <= 0.0
        return
    end
    for cluster=1:n_clusters
        if cluster_sizes[ cluster ] == 1
            clustering_coefficients[ cluster ] = 1.0
        else
            for i=2:cluster_sizes[ cluster ] - 1
                for j=i+1:cluster_sizes[cluster]
                    dist=distance_matrix[ assignments[ cluster, i], assignments[cluster,j] ]
                    if  dist > 0 && dist < cut_off
                        clustering_coefficients += 1
                    end
                end
            end
        end
        clustering_coefficients[ cluster ] /= cluster_sizes[i]*(cluster_sizes[i]-1)/2
    end
    return clustering_coefficients
end



# Daura Clustering
#==============================================================================#
# algorithm from Daura et al, Angew. Chem. Int. Ed. 38, 236-240, 1999
function dauraClustering( distance_matrix::Array{T2,2} , cut_off::T3 ) where { T1 <: Int, T2 <: Real, T3 <: Real }
	n_elements=size(distance_matrix)[1]
    cluster_centers=zeros(Int,n_elements)
	index_data=zeros(Int,n_elements)

	# Used
	used=zeros(Int,n_elements)  # vector to
	n_element_left = n_elements # number of elements left to assign
    n_clusters = 0              # number of clusters

	while sum(used) < n_elements
		# Basic info
		cluster_center=0
		nb_neighbor_max=0
		n_clusters += 1

		# Looking up
		for i=1:n_elements
			n_neighbor=0
			if used[i] != 0
				continue
			end
			for j=1:n_elements
				if used[j] != 0 || i == j
					continue
				end
				if distance_matrix[i,j] < cut_off
					n_neighbor +=1
				end
			end
			if nb_neighbor_max < n_neighbor
				nb_neighbor_max = n_neighbor
				cluster_center = i
			end
		end

		if cluster_center != 0
			# Assign clusters
			used[cluster_center] = -n_clusters
			index_data[cluster_center] = n_clusters
			for i=1:n_elements
				if i != cluster_center && used[i] == 0
					if distance_matrix[cluster_center,i] < cut_off
						used[i] = 1
						index_data[i] = n_clusters
					end
				end
			end
		else
			n_clusters -= 1
			break
		end
    end

	# Cluster centers and cluster sizes
	cluster_sizes=zeros(Int,n_clusters)
	cluster_centers=zeros(Int,n_clusters)
	count_cluster=1
	for i=1:n_elements
		if used[i] < 0
			cluster_nb=Int(abs(used[i]))
			cluster_centers[count_cluster] = i
			for j=1:n_elements
				if index_data[j]  == cluster_nb
					cluster_sizes[count_cluster] += 1
				end
			end
			count_cluster += 1
		end
	end

    # Return
    return cluster_centers, cluster_sizes, index_data
end
#==============================================================================#

# Kernels
#==============================================================================#
function gaussianKernel( matrix_distance::Array{T1,2} , cut_off_distance::T2) where { T1 <: Real, T2 <: Real }
    nb_element=size(matrix_distance)[1]
    rho=zeros(nb_element)
    for i=1:nb_element
        for j=i+1:nb_element
            gauss = exp(-(matrix_distance[i,j]/cut_off_distance)*(matrix_distance[i,j]/cut_off_distance))
            rho[i] += gauss
            rho[j] += gauss
        end
    end
    return rho
end
#==============================================================================#

# MAX
#==============================================================================#
function maxArray( matrix::Array{T1,2} ) where { T1 <: Real }
    nb_element=size(matrix)[1]
    max=0
    for i=1:nb_element
        for j=i+1:nb_element
            if matrix[i,j] > max
                max=matrix[i,j]
            end
        end
    end
    return max
end
function maxVector( vector::Vector{T1} ) where { T1 <: Real }
    max=0
    for i=1:size(vector)[1]
        if max < vector[i]
            max=vector[i]
        end
    end
    return max
end
#==============================================================================#

# Sequence vector
#==============================================================================#
function simpleSequence( size_::T1 ) where { T1 <: Int }
    vector=zeros(Int, size_)
    for i=1:size_
        vector[i] = i
    end
    return vector
end
#==============================================================================#

# Sort
#==============================================================================#
function sortMatrix( x::Vector{T1} ) where { T1 <: Real }
    sizex=size(x)[1]
    index_x=zeros(sizex)
    for i=1:sizex
        index_x[i]=i
    end
    for i=1:sizex
        for j=i:sizex
            if x[i] > x[j]
                stock=x[i]
                stocki=index_x[i]
                x[i]=x[j]
                index_x[i]=index_x[j]
                x[j]=stock
                index_x[j]=stocki
            end
        end
    end
    return index_x
end
function sortDescend( vector::Vector{T1} ) where { T1 <: Real }
    size_vector=size(vector)[1]
    vector_sorted=copy(vector)
    index_vector=simpleSequence(size_vector)
    for i=1:size_vector
        for j=i+1:size_vector
            if vector_sorted[i] < vector_sorted[j]
                stock=vector_sorted[i]
                stock2=index_vector[i]
                vector[i]=vector_sorted[j]
                index_vector[i]=index_vector[j]
                vector_sorted[j]=stock
                index_vector[j]=stock2
            end
        end
    end
    return vector_sorted, index_vector
end
function sortDescend( data::Array{T1,2}, index::T2 )  where { T1 <: Real , T2 <: Int }
	size_data=size(data)[1]
	n_dim=size(data)[2]
	for i=1:size_data-1
		for j=i+1:size_data
			if data[i,index] < data[j,index]
				stock = data[i,:]
				data[i,:] = data[j,:]
				data[j,:] = stock
			end
		end
	end
	return
end
#==============================================================================#

# Density Peak Clustering
#==============================================================================#
# Algorithm from Rodriguez, Alex, and Alessandro Laio. "Clustering by fast search and find of density peaks." Science 344.6191 (2014): 1492-1496.
function densityPeakClusteringTrain( data::Array{T1,2}, dc::T2 ) where { T1 <: Real, T2 <: Real }

	# Size and dimension of the input data
	size_data = size(data)[1]
	n_dim = size(data)[2]
	min_delta=0.1  # Decision min-delta to be cluster center
	min_rho=0.1    # Decision min-rho to be cluster center

	# Compute the maximum values of each dimensions
	max_v=data[1,:]
	min_v=data[1,:]
	for i=1:size_data
	    for j=1:n_dim
	        if max_v[j] < data[i,j]
	            max_v[j] = data[i,j]
	        end
	        if min_v[j] > data[i,j]
	            min_v[j] = data[i,j]
	        end
	    end
	end

	# Compute the distance matrix - most computation expensive
	# and memory consuming part; each dimension is normalized between 0 and 1
	distance_matrix, max_distance = computeDistanceMatrixAndMax( data , n_dim, max_v, min_v )

	# Compute the rho (~ local density of each points)
	rho = gaussianKernel( distance_matrix, dc)

	# Sort the rhos in decreasing order
	index=simpleSequence(size(rho)[1])
	for i=1:size(rho)[1]
		for j=i+1:size(rho)[1]
			if rho[i] < rho[j]
				stock=rho[i]
				rho[i]=rho[j]
				rho[j]=stock
				stock=index[i]
				index[i]=index[j]
				index[j]=stock
			end
		end
	end

	# Compute delta
	delta=ones(size_data)*max_distance
	delta[ 1 ] = -1
	nneigh=zeros(Int,size_data)
	for i=2:size_data
		for j=1:i-1
			if distance_matrix[ index[i] , index[j] ] < delta[ i ]
				delta[ i ] = distance_matrix[ index[i], index[j] ]
				nneigh[ i ] = j
			end
		end
	end

	# Compute maximum rho
	max_rho=rho[1]

	# put the delta of the max rho point to the max of delta
	for i=1:size_data
		if distance_matrix[index[1],i] > delta[1]
			delta[1] = distance_matrix[index[1],i]
		end
	end

	# Computing the max delta and assigning it to the largest cluster centers
	max_delta=0
	for i=1:size(delta)[1]
		if delta[i] > max_delta
			max_delta  = delta[i]
		end
	end

	n_cluster=0 # Number of clusters
	icl=[]      # Index of cluster centers
	cl=ones(Int,size_data)*(-1) # Assignement of the data points (cluster #)

	# Determine the cluster centers
	for i=1:size_data
	    if rho[i]/max_rho > min_rho && delta[i]/max_delta > min_delta
			n_cluster += 1
	        cl[index[i]] = n_cluster
	        icl=push!(icl,index[i])
	    end
	end

	if n_cluster != 0
	# Affectation of points to clusters using their nearest neighbor
		for i=1:size_data
	    	if cl[index[i]] == -1
	        	cl[index[i]] = cl[ index[nneigh[i]]  ]
	    	end
		end
		return cl, icl
	else
		return [], []
	end
end
function densityPeakClusteringTrain( data::Array{T1,2}, dc::T2 , min_rho::T3, min_delta::T4 ) where { T1 <: Real, T2 <: Real, T3 <: Real, T4 <: Real }

	# Size and dimension of the input data
	size_data = size(data)[1]
	n_dim = size(data)[2]
	min_delta=0.1  # Decision min-delta to be cluster center
	min_rho=0.1    # Decision min-rho to be cluster center

	# Compute the maximum values of each dimensions
	max_v=data[1,:]
	min_v=data[1,:]
	for i=1:size_data
	    for j=1:n_dim
	        if max_v[j] < data[i,j]
	            max_v[j] = data[i,j]
	        end
	        if min_v[j] > data[i,j]
	            min_v[j] = data[i,j]
	        end
	    end
	end

	# Compute the distance matrix - most computation expensive
	# and memory consuming part; each dimension is normalized between 0 and 1
	distance_matrix, max_distance = computeDistanceMatrixAndMax( data , n_dim, max_v, min_v )

	# Compute the rho (~ local density of each points)
	rho = gaussianKernel( distance_matrix, dc)

	# Sort the rhos in decreasing order
	index=simpleSequence(size(rho)[1])
	for i=1:size(rho)[1]
		for j=i+1:size(rho)[1]
			if rho[i] < rho[j]
				stock=rho[i]
				rho[i]=rho[j]
				rho[j]=stock
				stock=index[i]
				index[i]=index[j]
				index[j]=stock
			end
		end
	end

	# Compute delta
	delta=ones(size_data)*max_distance
	delta[ 1 ] = -1
	nneigh=zeros(Int,size_data)
	for i=2:size_data
		for j=1:i-1
			if distance_matrix[ index[i] , index[j] ] < delta[ i ]
				delta[ i ] = distance_matrix[ index[i], index[j] ]
				nneigh[ i ] = j
			end
		end
	end

	# Compute maximum rho
	max_rho=0
	for i=1:size(rho)[1]
		if rho[i] > max_rho
			max_rho = rho[i]
		end
	end

	# put the delta of the max rho point to the max of delta
	for i=1:size_data
		if distance_matrix[index[1],i] > delta[1]
			delta[1] = distance_matrix[index[1],i]
		end
	end

	# Computing the max delta and assigning it to the largest cluster centers
	max_delta=0
	for i=1:size(delta)[1]
		if delta[i] > max_delta
			max_delta  = delta[i]
		end
	end

	n_cluster=0 # Number of clusters
	icl=[]      # Index of cluster centers
	cl=ones(Int,size_data)*(-1) # Assignement of the data points (cluster #)

	# Determine the cluster centers
	for i=1:size_data
	    if rho[i] > min_rho && delta[i] > min_delta
			n_cluster += 1
	        cl[index[i]] = n_cluster
	        icl=push!(icl,index[i])
	    end
	end

	# Affectation of points to clusters using their nearest neighbor
	for i=1:size_data
	    if cl[index[i]] == -1
	        cl[index[i]] = cl[ index[nneigh[i]]  ]
	    end
	end

	return cl, icl
end
function densityPeakClusteringFirstStepDistanceMatrix( distance_matrix::Array{T1,2}, dc::T2 , file::T3 ) where { T1 <: Real, T2 <: Real, T3 <: AbstractString }

	size_data=size(distance_matrix)[1]

	max_distance=0
	for i=1:size_data-1
		for j=i+1:size_data
			if max_distance < distance_matrix[i,j]
				max_distance = distance_matrix[i,j]
			end
		end
	end

	# Compute the rho (~ local density of each points)
	rho = gaussianKernel( distance_matrix, dc)

	# Sort the rhos in decreasing order
	index=simpleSequence(size(rho)[1])
	for i=1:size(rho)[1]
		for j=i+1:size(rho)[1]
			if rho[i] < rho[j]
				stock=rho[i]
				rho[i]=rho[j]
				rho[j]=stock
				stock=index[i]
				index[i]=index[j]
				index[j]=stock
			end
		end
	end

	# Compute delta
	delta=ones(size_data)*max_distance
	delta[ 1 ] = -1
	nneigh=zeros(Int,size_data)
	for i=2:size_data
		for j=1:i-1
			if distance_matrix[ index[i] , index[j] ] < delta[ i ]
				delta[ i ] = distance_matrix[ index[i], index[j] ]
				nneigh[ i ] = j
			end
		end
	end

	# Compute maximum rho
	max_rho=0
	for i=1:size(rho)[1]
		if rho[i] > max_rho
			max_rho = rho[i]
		end
	end

	# put the delta of the max rho point to the max of delta
	for i=1:size_data
		if distance_matrix[index[1],i] > delta[1]
			delta[1] = distance_matrix[index[1],i]
		end
	end

	# Computing the max delta and assigning it to the largest cluster centers
	max_delta=0
	for i=1:size(delta)[1]
		if delta[i] > max_delta
			max_delta  = delta[i]
		end
	end

	# Writting decision diagram
	file_out=open(file,"w")
	for i=1:size(rho)[1]
		write(file_out,string(rho[i]," ",delta[i],"\n"))
	end
	close(file_out)

	return rho, delta, index, nneigh
end
function densityPeakClusteringFirstStep( data::Array{T1,2}, dc::T2 , file::T3 ) where { T1 <: Real, T2 <: Real, T3 <: AbstractString }

	# Size and dimension of the input data
	size_data = size(data)[1]
	n_dim = size(data)[2]
	min_delta=0.1  # Decision min-delta to be cluster center
	min_rho=0.1    # Decision min-rho to be cluster center

	# Compute the maximum values of each dimensions
	max_v=data[1,:]
	min_v=data[1,:]
	for i=1:size_data
	    for j=1:n_dim
	        if max_v[j] < data[i,j]
	            max_v[j] = data[i,j]
	        end
	        if min_v[j] > data[i,j]
	            min_v[j] = data[i,j]
	        end
	    end
	end

	# Compute the distance matrix - most computation expensive
	# and memory consuming part; each dimension is normalized between 0 and 1
	distance_matrix, max_distance = computeDistanceMatrixAndMax( data , n_dim, max_v, min_v )

	# Compute the rho (~ local density of each points)
	rho = gaussianKernel( distance_matrix, dc)

	# Sort the rhos in decreasing order
	index=simpleSequence(size(rho)[1])
	for i=1:size(rho)[1]
		for j=i+1:size(rho)[1]
			if rho[i] < rho[j]
				stock=rho[i]
				rho[i]=rho[j]
				rho[j]=stock
				stock=index[i]
				index[i]=index[j]
				index[j]=stock
			end
		end
	end

	# Compute delta
	delta=ones(size_data)*max_distance
	delta[ 1 ] = -1
	nneigh=zeros(Int,size_data)
	for i=2:size_data
		for j=1:i-1
			if distance_matrix[ index[i] , index[j] ] < delta[ i ]
				delta[ i ] = distance_matrix[ index[i], index[j] ]
				nneigh[ i ] = j
			end
		end
	end

	# Compute maximum rho
	max_rho=0
	for i=1:size(rho)[1]
		if rho[i] > max_rho
			max_rho = rho[i]
		end
	end

	# put the delta of the max rho point to the max of delta
	for i=1:size_data
		if distance_matrix[index[1],i] > delta[1]
			delta[1] = distance_matrix[index[1],i]
		end
	end

	# Computing the max delta and assigning it to the largest cluster centers
	max_delta=0
	for i=1:size(delta)[1]
		if delta[i] > max_delta
			max_delta  = delta[i]
		end
	end

	# Writting decision diagram
	file_out=open(file,"w")
	for i=1:size(rho)[1]
		write(file_out,string(rho[i]/max_rho," ",delta[i]/max_delta,"\n"))
	end
	close(file_out)

	return rho, delta
end
function densityPeakClusteringFirstStep( data::Array{T1,2}, dc::T2 , min_rho::T3, min_delta::T4, file::T5 ) where { T1 <: Real, T2 <: Real, T3 <: Real, T4 <: Real, T5 <: AbstractString }

	# Size and dimension of the input data
	size_data = size(data)[1]
	n_dim = size(data)[2]
	min_delta=0.1  # Decision min-delta to be cluster center
	min_rho=0.1    # Decision min-rho to be cluster center

	# Compute the maximum values of each dimensions
	max_v=data[1,:]
	min_v=data[1,:]
	for i=1:size_data
	    for j=1:n_dim
	        if max_v[j] < data[i,j]
	            max_v[j] = data[i,j]
	        end
	        if min_v[j] > data[i,j]
	            min_v[j] = data[i,j]
	        end
	    end
	end

	# Compute the distance matrix - most computation expensive
	# and memory consuming part; each dimension is normalized between 0 and 1
	distance_matrix, max_distance = computeDistanceMatrixAndMax( data , n_dim, max_v, min_v )

	# Compute the rho (~ local density of each points)
	rho = gaussianKernel( distance_matrix, dc)

	# Sort the rhos in decreasing order
	index=simpleSequence(size(rho)[1])
	for i=1:size(rho)[1]
		for j=i+1:size(rho)[1]
			if rho[i] < rho[j]
				stock=rho[i]
				rho[i]=rho[j]
				rho[j]=stock
				stock=index[i]
				index[i]=index[j]
				index[j]=stock
			end
		end
	end

	# Compute delta
	delta=ones(size_data)*max_distance
	delta[ 1 ] = -1
	nneigh=zeros(Int,size_data)
	for i=2:size_data
		for j=1:i-1
			if distance_matrix[ index[i] , index[j] ] < delta[ i ]
				delta[ i ] = distance_matrix[ index[i], index[j] ]
				nneigh[ i ] = j
			end
		end
	end

	# Compute maximum rho
	max_rho=0
	for i=1:size(rho)[1]
		if rho[i] > max_rho
			max_rho = rho[i]
		end
	end

	# put the delta of the max rho point to the max of delta
	for i=1:size_data
		if distance_matrix[index[1],i] > delta[1]
			delta[1] = distance_matrix[index[1],i]
		end
	end

	# Computing the max delta and assigning it to the largest cluster centers
	max_delta=0
	for i=1:size(delta)[1]
		if delta[i] > max_delta
			max_delta  = delta[i]
		end
	end

	# Writting decision diagram
	file_out=open(file,"w")
	for i=1:size(rho)[1]
		write(file_out,string(rho[i]," ",delta[i],"\n"))
	end
	close(file_out)

	return rho, delta
end
function densityPeakClusteringSecondStep( rho::Vector{T1}, delta::Vector{T2}, index::Vector{T3}, nearest_neighbor::Vector{T4}, min_rho::T5, min_delta::T6 ) where { T1 <: Real, T2 <: Real, T3 <: Int, T4 <: Int, T5 <: Real, T6 <: Real }

	size_data=size(rho)[1]
	n_cluster=0 # Number of clusters
	icl=[]      # Index of cluster centers

	cl=ones(Int,size_data)*(-1) # Assignement of the data points (cluster #)

	# Determine the cluster centers
	for i=1:size_data
	    if rho[i] > min_rho && delta[i] > min_delta
			n_cluster += 1
	        cl[index[i]] = n_cluster
	        icl=push!(icl,index[i])
	    end
	end

	# Affectation of points to clusters using their nearest neighbor
	for i=1:size_data
	    if cl[index[i]] == -1
	        cl[index[i]] = cl[ index[nearest_neighbor[i]]  ]
	    end
	end

	return cl, icl
end
#==============================================================================#

end
