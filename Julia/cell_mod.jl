module cell_mod

export Cell_param, Cell_vec, Cell_matrix
export Cell
export vec2matrix, wrap, dist1D, distance, compressParams, compressAtoms
export velocityFromPosition

# Import all import module
#----------------------------
using LinearAlgebra
using atom_mod
#----------------------------

#-------------
# Structures
#-----------------------------
mutable struct Cell_param
    length::Vector{Real}
    angles::Vector{Real}
    function Cell_param()
        new(zeros(Real,3),zeros(Real,3));
    end
    function Cell_param( params::Vector{T1} ) where { T1 <: Real }
        new(params,ones(Real,3)*90)
    end
    function Cell_param( a::T1, b::T2, c::T3 ) where { T1 <: Real, T2<: Real, T3 <: Real }
        new([a,b,c],[90.,90.,90.])
    end
    function Cell_param( a::T1, b::T2, c::T3, alpha::T4, beta::T5, gamma::T6 ) where { T1 <: Real, T2<: Real, T3 <: Real, T4 <: Real, T5 <: Real, T6 <: Real }
        new([a,b,c],[alpha,beta,gamma])
    end
    function Cell_param( lengths::Vector{T1}, angles::Vector{T2} ) where { T1 <: Real, T2 <: Real }
        new( lengths, angles )
    end
end
# Probably will be deleted soon as it is useless
mutable struct Cell_vec
    v1::Vector{Real}
    v2::Vector{Real}
    v3::Vector{Real}
    function Cell_vec( )
        new([],[],[]);
    end
    function Cell_vec( v1::Vector{T1} ) where { T1 <: Real }
        if size(v1)[1] == 3
            v2    = zeros(3)
            v3    = zeros(3)
            v2[2] = v1[1]
            v3[3] = v1[1]
            new( v1, v2, v3)
        else
            print("Error building vector !\n")
            new( [], [], [] )
        end
    end
    function Cell_vec( v1::Vector{T1}, v2::Vector{T2}, v3::Vector{T3} ) where { T1 <: Real, T2 <: Real, T3 <: Real }
        if size(v1)[1] == 3 && size(v2)[1] == 3 && size(v3)[1] == 3
            new( v1, v2, v3 )
        else
            print("Error building vector !\n")
            new( [], [], [] )
        end
    end
end
mutable struct Cell_matrix
    matrix::Array{Real,2}
    function Cell_matrix()
        new(Array{Real}(3,3));
    end
    function Cell_matrix( matrix::Array{T1}) where { T1 <: Real }
        if size(matrix)[1]==3 && size(matrix)[2]==3
            new( matrix )
        end
    end
    function Cell_matrix( a::T1, b::T1, c::T1 ) where { T1 <: Real }
        new( [ [ a, 0, 0 ], [ 0, b, 0 ] , [ 0, 0, c ] ] )
    end
    function Cell_matrix( v1::Vector{T1}, v2::Vector{T2}, v3::Vector{T3} ) where { T1 <: Real, T2 <: Real, T3 <: Real }
        if size(v1)[1] == 3 && size(v2)[1] == 3 && size(v3)[1] == 3
            new( [ v1, v2, v3 ] )
        else
            print("Error building matrix!\n")
            new( Array{Real}( 3, 3 ) )
        end
    end
end
#------------------------------

#------------------------------
# General type and conversions
#--------------------------------------------
Cell=Union{Cell_param, Cell_vec, Cell_matrix}
#---------------------------------------------

# Conversions
#-------------------------------------------------------------------------------
function cellMatrix2Params( cell_matrix::Array{T1,2} )  where { T1 <: Real }
    params=zeros(Real,3)
    for i=1:size(cell_matrix)[1]
        for j=1:size(cell_matrix)[2]
            params[i] += cell_matrix[i,j]*cell_matrix[i,j]
        end
        params[i] = sqrt( params[i])
    end
    return params
end
function cellMatrix2Params( cell_matrix::T1 )  where { T1 <: Cell_matrix }
    return cellMatrix2Params(cell_matrix.matrix)
end
function cellVector2Matrix( vectors::T1 ) where { T1 <: Cell_vec }
    matrix=Cell_matrix()
    for i=1:3
        matrix.matrix[i,1] = vectors.v1[i]
    end
    for i=1:3
        matrix.matrix[i,2] = vectors.v2[i]
    end
    for i=1:3
        matrix.matrix[i,3] = vectors.v3[i]
    end
    return matrix
end
function params2Matrix( cell_params::T1 ) where { T1 <: Cell_param }
    matrix=zeros(3,3)
    lengths=copy(cell_params.length)
    for i=2:3
        lengths[i] /= lengths[1]
    end
    tau=pi/180
    # COS of angles
    cos_ang=zeros(3)
    for i=1:3
        cos_ang[i] = cos( cell_params.angles[i]*tau)
    end
    temp=sqrt(1.0-cos_ang[3]*cos_ang[3])
    matrix[1,1] = lengths[1]
    matrix[1,2] = lengths[1]*lengths[2]*cos_ang[3]
    matrix[2,2] = lengths[1]*lengths[2]*temp
    matrix[3,1] = lengths[1]*lengths[3]*cos_ang[2]
    matrix[3,2] = lengths[1]*lengths[3]*(cos_ang[1]-cos_ang[2]*cos_ang[3])/temp

    temp=(1.0 + 2.0 *cos_ang[1]*cos_ang[2]*cos_ang[3] - cos_ang[1]*cos_ang[1] - cos_ang[2]*cos_ang[2] - cos_ang[3]*cos_ang[3])
    matrix[3,3] = lengths[1]*lengths[3]*sqrt(temp/(1-cos_ang[3]*cos_ang[3]))
    return Cell_matrix(matrix)
end
#---------------------------------------------------------------------------\

#-------------------------------------------------------------------------------
function getVolume( cell_matrix::Array{T1,2}) where { T1 <: Real }
    return LinearAlgebra.det(cell_matrix)
end
function getVolume( cell_matrix::T1) where { T1 <: Cell_matrix }
    return LinearAlgebra.det(cell_matrix.matrix)
end
function getVolume( cell_param::T1) where { T1 <: Cell_param }
    return LinearAlgebra.det(params2Matrix(cell_param))
end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
function wrap( position::T1, length::T2 ) where { T1 <: Real, T2 <: Real}
    sign=-1
    if position < 0
        sign=1
    end
    while position < 0 || position > length
        position = position + sign*length
    end
    return position
end
function wrap( atoms::T1, cell::T2 ) where { T1 <: atom_mod.AtomList, T2 <: Cell_matrix }
    # Computes cell parameters
    #--------------------------------------------
    params=[0.,0.,0.]
    for i=1:3
        for j=1:3
            params[i]=params[i]+cell.matrix[i,j]
        end
    end
    #--------------------------------------------

    #---------------
    # Compute atoms
    #---------------------------------
    for i=1:size(atoms.positions)[1]
        for j=1:3
            atoms.positions[i,j] = wrap( atoms.positions[i,j],params[j])
        end
    end
    #----------------------------------

    return atoms
end
function wrap( atoms::T1, cell::T2 ) where { T1 <: atom_mod.AtomList, T2 <: Cell_param }
    for i=1:size(atoms.positions)[1]
        for j=1:3
            atoms.positions[i,j] = wrap( atoms.positions[i,j],cell.length[j])
        end
    end
    return atoms
end
function wrap( molecules::T1, cell::T2 ) where { T1 <: atom_mod.AtomMolList, T2 <: Cell_param }
    for i=1:size(molecules.positions)[1]
        for j=1:3
            molecules.positions[i,j] = wrap( molecules.positions[i,j],cell.length[j])
        end
    end
    return molecules
end
function wrap( positions::Array{T1,2}, cell::T2 ) where { T1 <: Real, T2 <: Cell_param }
    for j=1:size(positions)[1]
        for i=1:3
            positions[j,i] = wrap( positions[j,i],cell.length[i] )
        end
    end
    return positions
end
function wrap( positions::Vector{T1}, cell::T2 ) where { T1 <: Real, T2 <: Cell_param }
    for i=1:3
        positions[i] = wrap( positions[i],cell.length[i] )
    end
    return positions
end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
function invertCell( cell_matrix::Array{T1,2} ) where { T1 <: Real }
    if det(cell_matrix) != 0
        return LinearAlgebra.inv(cell_matrix)
    else
        return false
    end
end
function invertCell( cell_matrix::T1 ) where { T1 <: Cell_matrix }
    return invertCell(cell_matrix.matrix)
end
function invertCell( cell_params::T1 ) where { T1 <: Cell_param }
    return invertCell(params2Matrix(cell_matrix.matrix))
end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
function getTransformedPosition( target_vector::Vector{T1}, cell_matrix::Array{T2,2} ) where {T1 <: Real, T2 <: Real }
    vector=zeros(3)
    for i=1:3
        for j=1:3
            vector[i] += cell_matrix[i,j]*target_vector[j]
        end
    end
    return vector
end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
function getScaledPosition( vector::Vector{T1}, cell_matrix::Array{T2,2} ) where { T1 <: Real , T2 <: Real }
    cell_mat_inverse=invertCell(cell_matrix)
    return getTransformedPosition(vector,cell_mat_inverse)
end
function getScaledPosition( vector::Vector{T1}, cell_matrix::T2 ) where { T1 <: Real , T2 <: Cell_matrix }
    return getScaledPosition(vector,cell_matrix.matrix)
end
function getScalePosition( target_matrix::Array{T1,2}, cell_matrix::Array{T2,2} ) where { T1 <: Real, T2 <: Real }
    inv_cell_matrix=invertCell(cell_matrix)
    nb_atoms=size(target_matrix)[1]
    scaled_positions=copy(target_matrix)
    for i=1:nb_atoms
        scaled_positions[i,:]=getTransformedPosition(target_matrix[i,:],inv_cell_matrix)
    end
    return scaled_positions
end
function getScalePosition( target_matrix::Array{T1,2}, cell_matrix::T2 ) where { T1 <: Real, T2 <: Cell_matrix }
    return getScalePosition(target_matrix,cell_matrix.matrix)
end
#-------------------------------------------------------------------------------

# Distance related functions
#-------------------------------------------------------------------------------
function dist1D( x1::T1, x2::T2, a::T3 ) where { T1 <: Real, T2 <: Real, T3 <: Real }
    dx=x1-x2
    if dx > a*0.5
        dx -= a
    end
    if dx < -a*0.5
        dx += a
    end
    return dx*dx
end
function distance( v1::Vector{T1}, v2::Vector{T2}, cell::Vector{T3} ) where { T1 <: Real, T2 <: Real, T3 <: Real }
    dist=0
    for i=1:size(v1)[1]
        dist += dist1D(v1[i],v2[i],cell[i])
    end
    return sqrt(dist)
end
function distance( v1::Vector{T1}, v2::Vector{T2}, cell_matrix::Array{T3,2} ) where { T1 <: Real, T2 <: Real, T3 <: Real }
    v1_scaled=getScaledPosition(v1,cell_matrix)
    v2_scaled=getScaledPosition(v2,cell_matrix)
    ds=zeros(3)
    # Distance + Min Image Convention
    for i=1:3
        ds[i] = v1_scaled[i]-v2_scaled[i]
        ds[i] = ds[i] - round(ds[i])
    end
    # Descaling of distance vector
    ds=getTransformedPosition(ds,cell_matrix)
    return sqrt(dot(ds,ds))
end
function distance( v1::Vector{T1}, v2::Vector{T2}, cell_matrix::Array{T3,2}, inv_cell_matrix::Array{T4,2} ) where { T1 <: Real, T2 <: Real, T3 <: Real, T4 <: Real }
    v1_scaled=getTransformedPosition(v1,inv_cell_matrix)
    v2_scaled=getTransformedPosition(v2,inv_cell_matrix)
    ds=zeros(3)
    # Distance + Min Image Convention
    for i=1:3
        ds[i] = v1_scaled[i]-v2_scaled[i]
        ds[i] = ds[i] - round(ds[i])
    end
    # Descaling of distance vector
    ds=getTransformedPosition(ds,cell_matrix)
    return sqrt(dot(ds,ds))
end
function distanceScale( v1_scaled::Vector{T1}, v2_scaled::Vector{T2}, cell_matrix::Array{T3,2} ) where { T1 <: Real, T2 <: Real, T3 <: Real }
    ds=zeros(3)
    # Distance + Min Image Convention
    for i=1:3
        ds[i] = v1_scaled[i]-v2_scaled[i]
        ds[i] = ds[i] - round(ds[i])
    end
    # Descaling of distance vector
    ds=getTransformedPosition(ds,cell_matrix)
    return sqrt(dot(ds,ds))
end
function distanceScale( v1_scaled::Vector{T1}, v2_scaled::Vector{T2}, cell_matrix::T3 ) where { T1 <: Real, T2 <: Real, T3 <: Cell_matrix }
    return distanceScale( v1_scaled, v2_scaled, cell_matrix )
end
function distance( v1::Vector{T1}, v2::Vector{T2}, cell_matrix::T3 ) where { T1 <: Real, T2 <: Real, T3 <: Cell_matrix }
    return distance(v1,v2,cell_matrix.matrix)
end
function distance( v1::Vector{T1}, v2::Vector{T2}, cell_params::T3 ) where { T1 <: Real, T2 <: Real, T3 <: Cell_param }
    return distance(v1,v2,params2Matrix(cell_params))
end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
function distance( positions::Array{T1,2}, cell::T2, atom1::T3, atom2::T4 ) where { T1 <: Real,  T2 <: Cell_param, T3 <: Real, T4 <: Real }
    dist=0
    for i=1:3
        dist +=  dist1D(positions[atom1,i],positions[atom2,i],cell.length[i])
    end
    return sqrt(dist)
end
function distance( positions1::Vector{T1}, positions2::Vector{T2}, cell_length::Vector{T3} ) where { T1 <: Real, T2 <: Real, T3 <: Real }
    dist=0
    for i=1:3
        dist += dist1D(positions1[i],positions2[i],cell_length[i])
    end
    return sqrt(dist)
end
function distance( atoms::T1, cell::T2, index1::T3, index2::T4 ) where { T1 <: atom_mod.AtomList, T2 <: Cell_param , T3 <: Int, T4 <: Int }
    dis=0
    for i=1:3
        dis += dist1D( atoms.positions[index1,i],atoms.positions[index2,i], cell.length[i] )
    end
    return sqrt(dis)
end
function distance( molecules::T1, cell::T2, index1::T3, index2::T4 ) where { T1 <: atom_mod.AtomMolList, T2 <: Cell_param , T3 <: Int , T4 <: Int }
    dist=0
    for i=1:3
        dist += dist1D( molecules.positions[index1,i],molecules.positions[index2,i], cell.length[i] )
    end
    return sqrt(dist)
end
function distance( atoms::T1, cell::T2, index1::T3, index2::T3, wrap::T4 ) where { T1 <: atom_mod.AtomList, T2 <: Cell_param,  T3 <: Int, T4 <: Bool }
    if (  wrap )
        wrap(atoms,cell)
    end
    return distance(atoms,cell,index1,index2)
end
#---------------------------------------------------------------------------

#----------
# Compress
#---------------------------------------------------------------------------
function compressParams( cell::T1, fracs::Vector{T2} ) where { T1 <: Cell_param, T2 <: Real }
    for i=1:3
        cell.lengths[i] *= fracs[i]
    end
    return cell
end
#---------------------------------------------------------------------------
function compressAtoms( atoms::T1 , cell::T2, fracs::Vector{T3} ) where { T1 <: atom_mod.AtomList, T2 <: Cell_param, T3 <: Real }
    for i=1:size(atoms.names)[1]
        for j=1:3
            atoms.positions[i,j] *= fracs[j]
        end
    end
    return atoms
end
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
function velocityFromPosition( traj::Vector{T1}, dt::T2, dx::T3 ) where { T1 <: atom_mod.AtomList, T2 <: Real, T3 <: Real }
    nb_atoms=size(traj[1].names)[1]
    nb_step=size(traj)[1]
    velocities=zeros(nb_step-1,nb_atoms,3)
    for step=1:nb_step-1
        for atom=1:nb_atoms
            for i=1:3
                velocities[step,atom,i]=(traj[step].positions[atom,i]-traj[step+1].positions[atom,i])/dt*dx
            end
        end
    end
    return velocities
end
#---------------------------------------------------------------------------

# Unwrap atoms target with regard to atom origin, using the cell information
# about the PBC
# The modifications are effected on the position Array and nothing is returned
#---------------------------------------------------------------------------
function unWrapOrtho!( positions::Array{T1,2}, origin::T2, target::T3, cell::T4  ) where { T1 <: Real, T2 <: Int, T3 <: Int, T4 <: Cell_param }
    for i=1:3
        dist=(positions[origin,i]-positions[target,i])
        if dist > cell.length[i]*0.5
            positions[target,i] += cell.length[i]
        elseif dist < -cell.length[i]*0.5
            positions[target,i] -= cell.length[i]
        end
    end
    return
end
#---------------------------------------------------------------------------

# Unwrap atoms target with regard to atom origin, using the cell information
# about the PBC
# The modifications are effected on the position structure AtomList and nothing is returned
#---------------------------------------------------------------------------
function unWrapOrtho!( structure::T1, origin::T2, target::T3, cell::T4  ) where { T1 <: atom_mod.AtomList, T2 <: Int, T3 <: Int, T4 <: Cell_param }
    for i=1:3
        dist=(structure.positions[origin,i]-structure.positions[target,i])
        if dist > cell.length[i]*0.5
            structure.positions[target,i] += cell.length[i]
        elseif dist < -cell.length[i]*0.5
            structure.positions[target,i] -= cell.length[i]
        end
    end
    return
end
#---------------------------------------------------------------------------

# Return two bools:
# - Is the molecule an infinite chain?
# - Was the chain exploration went ok?
#---------------------------------------------------------------------------
function isInfiniteChain( visited::Vector{T1}, matrix::Array{T2,2}, adjacency_table::Vector{T3}, positions::Array{T4,2} , cell::T5, target::T6, index_atoms::Vector{T7}, cut_off::T8 ) where { T1 <: Int, T2 <: Real, T3 <: Any, T4 <: Real, T5 <: cell_mod.Cell_param, T6 <: Int, T7 <: Int, T8 <: Real }
    visited[target]=1
    nb_neighbor=size(adjacency_table[target])[1]
    for neigh=1:nb_neighbor
        if visited[adjacency_table[target][neigh]] == 0
            unWrapOrtho!( positions, index_atoms[target], index_atoms[ adjacency_table[target][neigh] ], cell )
            if geom.distance( positions[ index_atoms[target], : ], positions[ index_atoms[ adjacency_table[target][neigh] ] , : ] ) > cut_off
                return false, false
            end
            isinf, isok = isInfiniteChain(visited,matrix,adjacency_table,positions,cell,adjacency_table[target][neigh],index_atoms,cut_off)
            # If infinite molecule is spotted, we stop
            if ! isok
                return true, false
            end
            if ! isinf
                return true, true
            end
        elseif geom.distance(positions[index_atoms[target],:],positions[ index_atoms[adjacency_table[target][neigh]] ,: ] ) > cut_off
            # Spotted infinite loop; stops the search
            return true, true
        end
    end
    return false, true
end
#---------------------------------------------------------------------------

# Return two bools:
# - Is the molecule an infinite chain?
# - Was the chain exploration went ok?
# Unwraps the molecule, even if infinite for visualisation purposes
#---------------------------------------------------------------------------
function checkInfinityAndUnWrap( visited::Vector{T1}, matrix::Array{T2,2}, adjacency_table::Vector{T3}, positions::Array{T4,2} , cell::T5, target::T6, index_atoms::Vector{T7}, cut_off::T8 ) where { T1 <: Int, T2 <: Real, T3 <: Any, T4 <: Real, T5 <: cell_mod.Cell_param, T6 <: Int, T7 <: Int, T8 <: Real }
    visited[target]=1
    nb_neighbor=size(adjacency_table[target])[1]
    for neigh=1:nb_neighbor
        if visited[adjacency_table[target][neigh]] == 0
            unWrapOrtho!( positions, index_atoms[target], index_atoms[ adjacency_table[target][neigh] ], cell )
            checkInfinityAndUnWrap(visited,matrix,adjacency_table,positions,cell,adjacency_table[target][neigh],index_atoms,cut_off)
        end
    end
    return
end
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
function reduced2Cartesian( positions_reduced::Vector{T1}, cell_matrix::Array{T2,2} ) where { T1 <: Real, T2 <: Real }
    new_positions=zeros(Real,3)
    for i=1:3 # x,y,z
        for j=1:3 # 1,2,3
            new_positions[i] += positions_reduced[j]*cell_matrix[i,j]
        end
    end
    return positions_reduced
end
function reduced2Cartesian( positions_reduced::Array{T1,2}, cell_matrix::Array{T2,2} ) where { T1 <: Real, T2 <: Real }
    nb_atoms=size(positions_reduced)[1]
    for atom = 1:nb_atoms
        positions_reduced[atom,:] = reduced2Cartesian( positions_reduced[atom,:], cell_matrix )
    end
    return positions_reduced
end
function reduced2Cartesian( positions_reduced::Array{T1,2}, cell_matrix::T2 ) where { T1 <: Real, T2 <: Cell_matrix }
    return reduced2Cartesian( positions_reduced, cell_matrix.matrix)
end
#---------------------------------------------------------------------------

end
