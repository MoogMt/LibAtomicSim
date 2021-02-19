module graph

# Load necessary modules from LibAtomicSim
using geom
using utils

# Description
# - Functions to deal with graphs theory and graph exploration

# Exporting functions
export searchGroupMember, groupsFromMatrix, getGroupMember, getGroupMemberAll, getGroupsFromMatrix
export getSizeTrees
export getAdjacent2Vertex, getAllAdjacentVertex
export extractMatrixForTree, extractAllMatrixForTrees


# Computing Groups
#--------------------------------------------------------------------------------
# Search all group members starting froma  given index
function searchGroupMember( matrix::Array{T1,2}, list::Vector{T2}, index::T3 , group_nb::T4 ) where { T1 <: Real , T2 <: Real , T3 <: Int , T4 <: Int }
    # Argument
    # - matrix: adjacency matrix describing the links in the graphs
    # - list: vector of int, containing group numbers for each edge
    # - index: index of target edge of the search
    # - group_nb: number of already found groups
    # Output
    # - list: updated list of the elements in the groupss

    # Compute number of edges in the adjacency matrix
    nb_vertex = size(matrix)[1]

    # Loop over edges
    for vertex=1:nb_vertex
        # Check that the edge is linked with the target edge
        if matrix[index,vertex] > 0
            # Check that the edge is not already affected
            # NB: there might be something better to do here with an else
            if list[vertex] == 0
                # If not affected, affects it to local group
                list[vertex] = group_nb
                # Launches a new local search
                list = searchGroupMember( matrix, list, vertex, group_nb )
            end
        end
    end

    # Returns the list of index belonging to each group
    return list
end
# Computes which connected graphs from adjacency matrix
function groupsFromMatrix( matrix::Array{T1,2} ) where { T1 <: Real }
    # Argument
    # - matrix: Adjacency matrix that contains all bonds between each vertex
    # Output
    # - nb_tree: number of independant trees in the matrix
    # - vertex_index: index of group of each vertex

    # Initialize number of trees to 0
    nb_tree = 0

    # Computes number of vertices
    nb_vertex = size( matrix )[1]

    # Initialize group to 0
    vertex_index = zeros(Int, nb_vertex )

    # Loop over vertices
    for vertex=1:nb_vertex
        # Check that the vertex has not been already affected to a given group
        if vertex_index[vertex] == 0
            # Increments number of trees
            nb_tree += 1
            # Look over all group members of current vertex
            vertex_index = searchGroupMember( matrix, vertex_index, vertex, nb_tree )
        end
    end

    # Returns number of trees, and index of all vertices
    return nb_tree, vertex_index
end
# Get the group members of a given group
function getGroupMember( target_index::T1, vertex_index::Vector{T2} ) where { T1 <: Int, T2 <: Int }
    # Argument
    # - target_index: index of the target tree
    # - vertex_index: index tree of all the vertices
    # Output
    # - members: vector of int, containing index of all vertices belonging to the graph

    # Initialize members list of vertices
    members=zeros(Int,0)

    # Get number of vertices
    size_ = size(vertex_index)[1]

    # Loop over vertices
    for vertex=1:size_
        # Check if the index of the vertex group is the one of the target group
        if target_index == vertex_index[vertex]
            # If so, add vertex to list
            push!( members, vertex )
        end
    end

    # Returns members of the group as a vector of index
    return members
end
# Get the group members of all groups
function getGroupMemberAll( nb_tree::T1, vertex_index::Vector{T2} ) where { T1 <: Int, T2 <: Int }
    # Argument
    # - nb_tree: number of groups contained in the system
    # - vertex_index: contains the group number of each vertex
    # Output
    # members_all: list of list containing all members of each groups

    # Initialize list of list
    members_all=[]

    # Loop over trees
    for tree=1:nb_tree
        # Adds the list of members of the tree to the lsit
        push!( members_all, getGroupMember( tree, vertex_index ) )
    end

    # Returns the list of list with all members of each groups
    return members_all
end
# Compute all group members starting from matrix
function getGroupsFromMatrix( matrix::Array{T1,2} ) where { T1 <: Real }
    # Argument
    # - matrix: adjacency matrix of the graph
    # Output:
    # - list of list containing members of each groups

    # Construct trees, get number of trees and the tree index (indicating the
    # group it belongs to) of each vertex
    nb_tree, vertex_index = groupsFromMatrix( matrix )

    # Returns list of list containing all vertices of each group
    return getGroupMemberAll( nb_tree, vertex_index )
end
#--------------------------------------------------------------------------------

# Computing size of trees
#--------------------------------------------------------------------------------
function getSizeTrees( vertex_index::Vector{T1} ) where { T1 <: Int }
    # Argument
    # - vertex_index: index of group of each vertex
    # Output
    # - sizes: vector, int, contains the sizes of the trees
    # - max: maximum size of the trees

    # Initialize tree_index as the
    tree_index = unique( vertex_index )

    # Get the number of trees
    nb_tree = size( tree_index )[1]

    # Initialize sizes as 0
    sizes = zeros(Int, nb_tree )

    # Get number of vertices
    nb_vertex = size( vertex_index )[1]

    # Loop over trees
    for tree=1:nb_tree
        # Loop over vertices
        for vertex=1:nb_vertex
            # If the vertex is that of the current tree, increments size counter
            if tree_index[tree] == vertex_index[vertex]
                sizes[tree] += 1
            end
        end
    end

    # Returns vector with sizes and maximum of sizes
    return sizes, maximum(sizes)
end
#--------------------------------------------------------------------------------

# Get Adjacency for Vertices
#--------------------------------------------------------------------------------
# Get vertices adjacent to a given vertex
function getAdjacent2Vertex( index::T1, matrix::Array{T2,2} ) where { T1 <: Int, T2 <: Real }
    # Argument
    # - index: index of the target index (int)
    # - matrix: adjacency matrix of the graph/system
    # Output
    # - adjacent_vertex: vector int, contains the index of the adjacent
    # vector to the target vertex

    # Initialize vector containing index of adjacency vertices
    adjacent_vertex=zeros(Int,0)

    # Get the number of vector
    nb_vertex=size(matrix)[1]

    # Loop over vertices
    for vertex=1:nb_vertex
        # If the vertex is neighbor to the target vertex, add the index of the
        # vertex to the vector
        if matrix[index,vertex] == 1
            push!( adjacent_vertex, vertex )
        end
    end

    # Returns the vector of index of ajdacent int
    return adjacent_vertex
end
# Get vertices adjacent to all vertices
function getAllAdjacentVertex( matrix::Array{T1,2} ) where { T1 <: Real }
    # Argument
    # - matrix: adjacency matrix of the graph/system
    # Output
    # - adjacent_vertex: list of adjacency vertex

    # Initialize list of adjacency vertices
    adjacent_vertex = []

    # Get number of vertices
    nb_vertex = size(matrix)[1]

    #  Loop over vertex
    for vertex=1:nb_vertex
        # Add list of the adjacency vertex to the list
        push!( adjacent_vertex, getAdjacent2Vertex( vertex, matrix ) )
    end

    # Returns list of adjacency vertex
    return adjacent_vertex
end
#--------------------------------------------------------------------------------

# Extract Matrix from Trees
#--------------------------------------------------------------------------------
# Extract tree matrix from general matrix
function extractMatrixForTree( matrix::Array{T1,2}, tree::Vector{T2} ) where { T1 <: Real, T2 <: Int }
    # Argument
    # - matrix: adjacency matrix of the whole system
    # - tree: target tree, with vectors of all index of the vertices in the tree
    # Output
    # - matrix_local: adjacency matrix of the tree

    # Get the size of the target tree
    size_tree=size(tree)[1]

    # Initialize tree matrix
    matrix_local=zeros(Int,size_tree,size_tree)

    # Loop over vertices of the tree 1
    for vertex1=1:size_tree
        # Loop over vertices of the tree 2
        for vertex2=vertex1+1:size_tree
            # Get the local value of the matrix for vertex1-vertex2 adjacency value
            matrix_local[vertex1,vertex2] = matrix[ tree[vertex1], tree[vertex2] ]

            # Use the fact that the matrix is symmetric
            matrix_local[vertex2,vertex1] = matrix_local[vertex1,vertex2]
        end
    end

    # Returns local matrix of the tree
    return matrix_local
end
# Extract all matrices for all trees
function extractAllMatrixForTrees( matrix::Array{T1,2}, trees::Vector{T2} ) where { T1 <: Real, T2 <: Any }
    # Argument
    # - matrix: general adjacency matrix of the system
    # - trees: list of vectors containing index of members of trees
    # Output
    # - matrices: matrices

    # Initialize list of matrices
    matrices=[]

    # Get the number of trees
    nb_tree=size(trees)[1]

    # Loop over trees
    for tree=1:nb_tree
        # Computes local matrix and adds it to the list
        push!( matrices, extractMatrixForTree(matrix, trees[tree] ) )
    end

    # list of local matrices for all trees
    return matrices
end
#--------------------------------------------------------------------------------
end
