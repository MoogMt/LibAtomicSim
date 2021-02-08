
export Cell_param, Cell_matrix
export Cell
export vec2matrix, wrap, dist1D, distance, compressParams, compressAtoms
export velocityFromPosition

# Import all import module
#----------------------------
using LinearAlgebra
using atom_mod
using geom
using graph
#----------------------------

#-----------------------------
# Cell_param : describes the cell, using lengths and angles
mutable struct Cell_param

    # Variables
    #-------------------------
    length::Vector{Real} # Lengths of the cell (angstroms)
    angles::Vector{Real} # Angles of the cell (degrees)
    #-------------------------

    # Constructors
    #--------------------------------------------------
    # Creates a default_cell params (default lengths=1A, default angles=90°)
    function Cell_param()
        # Construct Object
        new( ones(Real,3), ones(Real,3)*90.0 );
    end
    # Create a orthorombic cell, based on lengths parameters
    function Cell_param( lengths::Vector{T1} ) where { T1 <: Real }
        # Arguments:
        # lengths: lengths of the
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
    #--------------------------------------------------
end
# Cell : contains cell information both in parameters and matrix form
mutable struct Cell

    # Object variables
    #-------------------------------
    lengths::Vector{Real} # Contains the lengths of the cell (a,b,c) in angstroms
    angles::Vector{Real}  # Angles of the cell (alpha,beta,gamma) in degrees
    matrix::Array{Real,2} # Cell matrix
    #-------------------------------

    # Constructor functions
    #-----------------------------------------------------
    # Create a default cell (lengths=1,angles=90)
    function Cell()
        # Create Object
        new( ones(Real,3), ones(3)*90.0, Matrix{Real}(I,3,3)*1.0 )
    end
    # Create an orthorombic cell using provided lengths
    function Cell( a::T1, b::T2, c::T3 ) where { T1 <: Real, T2 <: Real, T3 <: Real }

        # Create matrix
        #------------------------
        matrix=zeros(3,3)
        matrix[1,1] = a
        matrix[2,2] = b
        matrix[3,3] = c
        #------------------------

        # Create Object
        new( [a,b,c], ones(3)*90.0, matrix )
    end
    # Create lengths using provided lengths and angles (in degrees)
    function Cell( a::T1, b::T2, c::T3, alpha::T4, beta::T5, gamma::T6 ) where { T1 <: Real, T2 <: Real, T3 <: Real, T4 <: Real, T5 <: Real, T6 <: Real }

        # Converts angles in radians
        #----------------------------
        alpha2 = alpha*pi/180
        beta2 = beta*pi/180
        gamma2 = gamma*pi/180
        #------------------------------

        # Creating the cell matrix
        #----------------------------------
        matrix = zeros(3,3)
        matrix[1,1] = a
        matrix[1,2] = b*cos( gamma2 )
        matrix[1,3] = c*cos( beta2 )
        matrix[2,2] = b*sin( gamma2 )
        matrix[2,3] = c*( cos( alpha2 ) - cos( beta2 )*cos( gamma2 ) )/sin( gamma2 )
        volume=sqrt( 1 + 2*cos(alpha2)*cos(beta2)*cos(gamma2) -cos(alpha2)^2 -cos(beta2)^2 -cos(gamma2)^2 )
        matrix[3,3] = c*volume/sin( gamma2 )
        #------------------------------------

        # Create matrix
        new( [a,b,c], [alpha,beta,gamma], matrix )
    end
    #----------------------------------------------------------
end
#--------------------------------------------------------------

# Conversions
#-------------------------------------------------------------------------------
function cellMatrix2Params( cell_matrix::Array{T1,2} )  where { T1 <: Real }
    length = zeros( Real, 3 )
    for col=1:3
        for line=1:3
            length[col] += cell_matrix[line,col]*cell_matrix[line,col]
        end
        length[col] = sqrt( length[col] )
    end
    tau=180/pi
    angles = zeros( Real, 3 )
    angles[1] = acos(sum( cell_matrix[:,2].*cell_matrix[:,3] )/(length[2]*length[3]))*tau
    angles[2] = acos(sum( cell_matrix[:,1].*cell_matrix[:,3] )/(length[1]*length[3]))*tau
    angles[3] = acos(sum( cell_matrix[:,1].*cell_matrix[:,2] )/(length[1]*length[2]))*tau
    return Cell_param( length, angles )
end
function cellMatrix2Params( cell_matrix::T1 )  where { T1 <: Cell_matrix }
    return cellMatrix2Params( cell_matrix.matrix )
end
function cellMatrix2Params( cell_matrices::Vector{T1} ) where { T1 <: Cell_matrix }
    nb_step = size( cell_matrices )[1]
    cells_params = Vector{ Cell_param }(undef, nb_step )
    for step=1:nb_step
        cells_params[step] = cellMatrix2Params( cell_matrices[step] )
    end
    return cells_params
end

function params2Matrix( cell_params::T1 ) where { T1 <: Cell_param }
    matrix=zeros(3,3)
    lengths=copy( cell_params.length)
    angles=copy(cell_params.angles)*pi/180
    matrix[1,1] = lengths[1]
    matrix[1,2] = lengths[2]*cos( angles[3] )
    matrix[1,3] = lengths[3]*cos( angles[2] )
    matrix[2,2] = lengths[2]*sin( angles[3] )
    matrix[2,3] = lengths[3]*( cos( angles[1] ) - cos( angles[2] )*cos( angles[3] ) )/sin( angles[3] )
    volume=sqrt( 1 + 2*cos(angles[1])*cos(angles[2])*cos(angles[3]) -cos(angles[1])^2 -cos(angles[2])^2 -cos(angles[3])^2 )
    matrix[3,3] = lengths[3]*volume/sin( angles[3] )
    return Cell_matrix( matrix )
end
function params2Matrix( cells_params::Vector{T1} ) where { T1 <: Cell_param }
    nb_step = size(cells_params)[1]
    cells_ = Vector{ Cell_matrix }( undef, nb_step )
    for step=1:nb_step
        cells_[step] = params2Matrix( cells_params[step] )
    end
    return cells_
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

function makeCellParams( lengths::Array{T1,2}, angles::Array{T2,2} ) where { T1 <: Real, T2 <: Real }
    nb_step = size(lengths)[1]
    cells = Vector{Cell_param}(undef, nb_step )
    for step=1:nb_step
        cells[step] = Cell_param( lengths[step,:], angles[step,:] )
    end
    return cells
end

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

    # Getting Scaled positions
    atoms.positions = getTransformedPosition( atoms.positions, LinearAlgebra.inv(cell.matrix ) )

    #---------------
    # Compute atoms
    #---------------------------------
    nb_atoms=size(atoms.names)[1]
    for atom=1:nb_atoms
        for i=1:3
            if atoms.positions[atom,i] < 1
                atoms.positions[atom,i] += 1
            end
            if atoms.positions[atom,i] > 1
                atoms.positions[atom,i] -= 1
            end
        end
    end
    #----------------------------------

    # Descaling
    atoms.positions = getTransformedPosition( atoms.positions, cell.matrix )

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
function wrap( traj::Vector{T1}, cell::T2 ) where { T1 <: atom_mod.AtomList, T2 <: Cell_param }
    nb_step=size(traj)[1]
    for step =1:nb_step
        traj[step] = wrap( traj[step], cell )
    end
    return traj
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
function getTransformedPosition( target_matrix::Array{T1,2}, cell_matrix::Array{T2,2} ) where {T1 <: Real, T2 <: Real }
    nb_point=size(target_matrix)[1]
    matrix_transformed = zeros( nb_point, 3 )
    for point=1:nb_point
        matrix_transformed[ point, : ] = getTransformedPosition( target_matrix[ point,: ], cell_matrix )
    end
    return matrix_transformed
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
# TODO Some bugs in the definition of distances, make sure everything is correct,
# when possible use longer names to define more precisely...
# TODO SORT THIS MESS
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
function distance( positions::Array{T1,2}, cell::T2, atom1::T3, atom2::T4 ) where { T1 <: Real,  T2 <: Cell_param, T3 <: Real, T4 <: Real }
    dist=0
    for i=1:3
        dist +=  dist1D(positions[atom1,i],positions[atom2,i],cell.length[i])
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

# Computes the velocities in a traj by finite elemnet method from
# the positions (least worst option)
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

#---------------------------------------------------------------------------------------------------------
# Unwraps a target molecule, even if infinite, unwraping atoms only once, exploring the molecule as a tree
# Useful for visualization
# Recursive.
function unWrapOrthoOnce( visited::Vector{T1}, matrix::Array{T2,2}, adjacency_table::Vector{T3}, positions::Array{T4,2} , cell::T5, target::T6, index_atoms::Vector{T7}) where { T1 <: Int, T2 <: Real, T3 <: Any, T4 <: Real, T5 <: Cell_param, T6 <: Int, T7 <: Int }
    visited[target]=1
    nb_neighbor=size(adjacency_table[target])[1]
    for neigh=1:nb_neighbor
        if visited[adjacency_table[target][neigh]] == 0
            unWrapOrtho!( positions, index_atoms[target], index_atoms[ adjacency_table[target][neigh] ], cell )
            unWrapOrthoOnce(visited,matrix,adjacency_table,positions,cell,adjacency_table[target][neigh],index_atoms)
        end
    end
    return
end
#---------------------------------------------------------------------------



#---------------------------------------------------------------------------
# Return two bools:
# - Is the molecule an infinite chain?
# - Was the chain exploration went ok?
function isInfiniteChain( visited::Vector{T1}, matrix::Array{T2,2}, adjacency_table::Vector{T3}, positions::Array{T4,2} , cell::T5, target::T6, index_atoms::Vector{T7}, cut_off::T8 ) where { T1 <: Int, T2 <: Real, T3 <: Any, T4 <: Real, T5 <: cell_mod.Cell_param, T6 <: Int, T7 <: Int, T8 <: Real }
    visited[target]=1
    nb_neighbor=size(adjacency_table[target])[1]
    for neigh=1:nb_neighbor
        if visited[adjacency_table[target][neigh]] == 0
            unWrapOrtho!( positions, index_atoms[target], index_atoms[ adjacency_table[target][neigh] ], cell )
            if geom.distance( positions[ index_atoms[target], : ], positions[ index_atoms[ adjacency_table[target][neigh] ] , : ] ) > cut_off
                return true, true
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
# Check if molecule is infinite
# by unwraping once all atoms following a graph exploration method
function checkInfiniteChain( matrix::Array{T1,2},  positions::Array{T2,2} , cell::T3, molecule_indexs::Vector{T4}, cut_off::T5 ) where { T1 <: Real, T2 <: Real, T3 <: Cell_param, T4 <: Int, T5 <: Real }
    adjacent_molecule = graph.getAllAdjacentVertex(matrix)
    size_molecule = size( matrix )[1]
    visited=zeros(Int,size_molecule)
    cell_mod.unWrapOrthoOnce( visited, matrix , adjacent_molecule, positions, cell, 1, molecule_indexs )
    for atom=1:size_molecule
        for adj=1:size(adjacent_molecule[atom])[1]
            if norm( positions[ molecule_indexs[atom], :] - positions[ molecule_indexs[adjacent_molecule[atom][adj]], :] ) > cut_off
                return true
            end
        end
    end
    # If we get there, the molecule isn't infinite
    return false
end
# Find the atoms that are bonded but whose
# bonds are cut by the one unwrap policty
function findUnlinked( matrix::Array{T1,2},  positions::Array{T2,2} , cell::T3, molecule_indexs::Vector{T4}, cut_off::T5  ) where { T1 <: Real, T2 <: Real, T3 <: Cell_param, T4 <: Int, T5 <: Real }
    list=Array{Int}(undef,0,2)
    size_molecule=size(matrix)[1]
    for atom=1:size_molecule
        for atom2=atom+1:size_molecule
            # matrix[i,j] correspond to the actual bond, norm(positions(i,:)-positions(j,:)) corresponds to the positions in unwrapped
            # so if one is 1 and the other not, we have an infinity loop
            if matrix[atom,atom2] == 1 && norm(positions[atom,:]-positions[atom2,:]) > cut_off
                if size(list)[1] == 0
                    list = vcat( list, [ molecule_indexs[atom] molecule_indexs[atom2] ] )
                else
                    check = true
                    # Checking that we haven't already added the pair
                    for i=1:size(list)[1]
                        if list[i,:] == [ molecule_indexs[atom], molecule_indexs[atom2] ]
                            check = false
                        end
                    end
                    # If ok we add
                    if check
                        list = vcat( list, [ molecule_indexs[atom] molecule_indexs[atom2] ] )
                    end
                end
            end
        end
    end
    # If we get thee, the molecule isn't infinite
    return list
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
    return new_positions
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

#-------------------------------------------------------------------------------
function computeMoveVector( index::Vector{T1}, cell_matrix::Array{T2,2} ) where { T1 <: Int, T2 <: Real }
    moveVector = zeros(Real,3)
    for i=1:3
        for j=1:3
            moveVector[i] += index[i]*cell_matrix[i,j]
        end
    end
    return moveVector
end
function growCell( cell::Array{T1,2}, n_grow::Vector{T2} ) where { T1 <: Real, T2 <: Int }
    cell2 = copy(cell)
    for i=1:3
        cell2[:,i] = cell[:,i]*n_grow[i]
    end
    return Cell_matrix(cell2)
end
function growCell( cell::T1 , n_grow::Vector{T2} ) where { T1 <: Cell_matrix, T2 <: Int }
    return growCell( cell.matrix, n_grow )
end
function growCell( cell::T1, n_grow::Vector{T2} ) where { T1 <: Cell_param, T2 <: Int }
    return growCell( params2Matrix(cell), n_grow )
end
function makeSuperCell( atoms::T1, cell_matrix::Array{T2,2}, n_grow::Vector{T3} ) where { T1 <: AtomList, T2 <: Real, T3 <: Int }
    nb_atoms_base = size(atoms.names)[1]
    nb_atoms_new = nb_atoms_base
    for i=1:3
        nb_atoms_new *= n_grow[i]
    end
    new_atoms = AtomList( nb_atoms_new )
    count_ = 1
    for dirx=0:n_grow[1]-1
        for diry=0:n_grow[2]-1
            for dirz=0:n_grow[3]-1
                move_box=[dirx,diry,dirz]
                # OLD IMPLEMENTATION
                moveVector = zeros(Real,3)
                for xyz=1:3
                    for v123=1:3
                        moveVector[xyz] += move_box[v123]*cell_matrix[xyz,v123]
                    end
                end
                for atom = 1:nb_atoms_base
                    new_atoms.positions[count_,:] = atoms.positions[atom,:] .+ moveVector[:]
                    new_atoms.names[count_] = atoms.names[atom]
                    new_atoms.index[count_] = 1
                    count_ += 1
                end
            end
        end
    end
    super_cell = growCell( cell_matrix, n_grow )
    atom_mod.sortAtomsByZ!(new_atoms)
    for i=1:nb_atoms_new
        new_atoms.index[i] = i
    end
    return new_atoms, cell_mod.cellMatrix2Params( super_cell )
end
function makeSuperCell!( traj::Vector{T1}, cell::Array{T2,2}, n_grow::Vector{T3}  ) where { T1 <: atom_mod.AtomList, T2 <: Real, T3 <: Int }
    nb_step = size(traj)
    for i=1:nb_step
        traj[step] = duplicateAtoms( traj[step], cell, n_grow )
    end
    return traj
end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
function nonOrtho2OrthoByCut( cell::T1 ) where { T1 <: cell_mod.Cell_matrix }
    lengths = zeros( Real, 3 )
    lengths[1] = LinearAlgebra.norm( cell.matrix[:,1] )
    for i=2:3
        for j=1:3
            lengths[j] += cell.matrix[i,j]
        end
    end
    return cell_mod.Cell_param( lengths )
end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
function makeCells( lengths::Array{T1,2}, angles::Array{T2,2} ) where { T1 <: Real, T2 <: Real }
    nb_step = size(lengths)[1]
    cells=Vector{ Cell_param }(undef, nb_step )
    for step=1:nb_step
        cells[step] = Cell_param( lengths[step,:], angles[step,:] )
    end
    return cells
end
#-------------------------------------------------------------------------------

end
