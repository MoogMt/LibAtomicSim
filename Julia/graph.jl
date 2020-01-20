module graph

using geom
using utils

export searchGroupMember, groupsFromMatrix, getSizeTree
export getGroupMember, getGroupMemberAll, getGroupsFromMatrix
export getAdjacent2Vertex, getAllAdjacentVertex 
export extractMatrixForTree, extractAllMatrixForTrees

# Computing Groups
#==============================================================================#
function searchGroupMember( matrix::Array{T1,2}, list::Vector{T2}, index::T3 , group_nb::T4 ) where { T1 <: Real , T2 <: Real , T3 <: Int , T4 <: Int }
    for i=1:size(matrix)[1]
        if matrix[index,i] > 0
            if list[i] == 0
                list[i]=group_nb
                list=searchGroupMember(matrix,list,i,group_nb)
            end
        end
    end
    return list
end
function groupsFromMatrix( matrix::Array{T1,2} ) where { T1 <: Real }
    nb_tree=0
    nb_vertex=size(matrix)[1]
    vertex_index=zeros(Int,nb_vertex)
    for i=1:nb_vertex
        if vertex_index[i] == 0
            nb_tree += 1
            vertex_index = searchGroupMember(matrix,vertex_index,i,nb_tree)
        end
    end
    return nb_tree, vertex_index
end
function getGroupMember( target_index::T1, vertex_index::Vector{T2} ) where { T1 <: Int, T2 <: Int }
    members=zeros(Int,0)
    size_ = size(vertex_index)[1]
    for i=1:size_
        if target_index == vertex_index[i]
            push!(members,i)
        end
    end
    return members
end
function getGroupMemberAll( nb_tree::T1, vertex_index::Vector{T2} ) where { T1 <: Int, T2 <: Int }
    members_all=[]
    for tree=1:nb_tree
        push!( members_all, getGroupMember(tree,vertex_index) )
    end
    return members_all
end
function getGroupsFromMatrix( matrix::Array{T1,2} ) where { T1 <: Real }
    nb_tree, vertex_index = groupsFromMatrix( matrix )
    return getGroupMemberAll( nb_tree, vertex_index )
end
#==============================================================================#

# Computing size of trees
#==============================================================================#
function getSizeTrees( vertex_index::Vector{T1} ) where { T1 <: Int }
    sizes=[]
    max=0
    tree_index=unique(vertex_index)
    for tree in tree_index
        size=0
        for i=1:size(vertex_index)[1]
            if tree == vertex_index
                size+=1
            end
        end
        if size > max
            max=size
        end
        push!(sizes,size)
    end
    return sizes, max
end
#==============================================================================#

# Get Adjacency for Vertices
#==============================================================================#
function getAdjacent2Vertex( index::T1, matrix::Array{T2,2} ) where { T1 <: Int, T2 <: Real }
    adjacent_vertex=zeros(Int,0)
    nb_vertex=size(matrix)[1]
    for i=1:nb_vertex
        if matrix[index,i] == 1
            push!(adjacent_vertex,i)
        end
    end
    return adjacent_vertex
end
function getAllAdjacentVertex( matrix::Array{T1,2} ) where { T1 <: Real }
    adjacent_vertex = []
    nb_matrix = size(matrix)[1]
    for i=1:nb_matrix
        push!( adjacent_vertex, getAdjacent2Vertex( i, matrix ) )
    end
    return adjacent_vertex
end
#==============================================================================#

# Extract Matrix from Trees
#==============================================================================#
function extractMatrixForTree( matrix::Array{T1,2}, tree::Vector{T2} ) where { T1 <: Real, T2 <: Int }
    size_tree=size(tree)[1]
    matrix_local=zeros(Int,size_tree,size_tree)
    for i=1:size_tree
        for j=i+1:size_tree
            matrix_local[i,j] = matrix[ tree[i], tree[j] ]
            matrix_local[j,i] = matrix_local[i,j]
        end
    end
    return matrix_local
end
function extractAllMatrixForTrees( matrix::Array{T1,2}, trees::Vector{T2} ) where { T1 <: Real, T2 <: Any }
    matrices=[]
    nb_tree=size(trees)[1]
    for tree=1:nb_tree
        push!( matrices, extractMatrixForTree(matrix, trees[tree] ) )
    end
    return matrices
end
#==============================================================================#

end
