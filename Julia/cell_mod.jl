module cell_mod

export Cell_param, Cell
export Cell
export wrap, dist1D, distance
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
        # Output:
        # Cell_param that describes the cell

        # Construct Object
        new( ones(Real,3), ones(Real,3)*90.0 );
    end
    # Create a orthorombic cell, based on lengths parameters given as a vector
    function Cell_param( lengths::Vector{T1} ) where { T1 <: Real }
        # Arguments:
        # - lengths: lengths of the cell in angstroms
        # Output:
        # - Cell_param object that describes the cell

        # Creating Object
        new( lengths, ones(Real,3)*90.0 )
    end
    # Create an orthorombic cell, based on lengths given as scalars
    function Cell_param( a::T1, b::T2, c::T3 ) where { T1 <: Real, T2 <: Real, T3 <: Real }
        # Arguments
        # - a,b,c: lengths of the cell (in angstroms)
        # Output:
        # - Cell_param object that describes the cell

        # Construct Cell_params
        new( [a,b,c], ones(3)*90.0 )
    end
    # Construct object using all parameters given as scalars
    function Cell_param( a::T1, b::T2, c::T3, alpha::T4, beta::T5, gamma::T6 ) where { T1 <: Real, T2<: Real, T3 <: Real, T4 <: Real, T5 <: Real, T6 <: Real }
        # Arguments
        # - a,b,c: lengths of the cell (in angstroms)
        # - alpha,beta, gamma: angles of the cells (in degrees)
        # Output
        # - Cell_param object that describes the cell

        # Constructing Object
        new( [a,b,c], [alpha,beta,gamma] )
    end
    # Construct Cell_params using all parameters given as vectors
    function Cell_param( lengths::Vector{T1}, angles::Vector{T2} ) where { T1 <: Real, T2 <: Real }
        # Arguments:
        # - lengths: vector contaning lengths (a,b,c) of the cell in angstroms
        # - angles: vector containg angles (alpha,beta,gamma) of the cell in degrees
        # Output
        # - a cell_param object describing the cell

        # Constructing Object
        new( lengths, angles )
    end
    # Construct a Cell_params trajectory from a set of matrix contaning lengths (angstroms) and angles (degrees) of the cells
    function Cell_param( lengths::Array{T1,2}, angles::Array{T2,2} ) where { T1 <: Real, T2 <: Real }
        # Arguments:
        # - lengths: vector containing the lengths of the cell (angstroms)
        # - angles: vector containing the angles of the cell (degrees)
        # Output:
        # - cell_params: vector containing cell_params for each step, describing the cell at time t

        # Get the number of step of the trajectory
        nb_step = size(lengths)[1]

        # Initialize output
        cells_params = Vector{ Cell_param }(undef, nb_step )

        # Loop over time
        for step=1:nb_step
            # Construct the cell parameters at time t
            cells_params[step] = Cell_param( lengths[step,:], angles[step,:] )
        end

        # Construct object
        return cells_params
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
        # Arguments:
        # a,b,c: lengths of the cell in angstroms

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
        # Arguments:
        # - a,b,c: lengths of the cell in angstroms
        # - alpha,beta,gamma: angles of the cell in degrees

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
    # NB: We probably need to add more functions, because the set that is
    # here is relatively limiting
    #----------------------------------------------------------
end
#--------------------------------------------------------------

# Matrix <-> Parameters conversion
#-------------------------------------------------------------------------------
# Converts a matrix into cell parameters (does not convert lengths)
function matrix2Params( cell_matrix::Array{T1,2} )  where { T1 <: Real }
    # Arguments:
    # cell-matrix: the cell matrix given as a ... matrix
    # Output: A Cell_Param object that corresponds to the cell describes by the matrix

    # Computing lengths (no conversion is attempted here)
    #---------------------------------------------------------------------
    length = zeros( Real, 3 )
    for col=1:3
        for line=1:3
            length[col] += cell_matrix[line,col]*cell_matrix[line,col]
        end
        length[col] = sqrt( length[col] )
    end
    #---------------------------------------------------------------------

    # Computing angles, converting to degrees
    #---------------------------------------------------------------------
    tau=180/pi
    angles = zeros(Real, 3 )
    angles[1] = acos( sum( cell_matrix[:,2].*cell_matrix[:,3] )/(length[2]*length[3]) )*tau
    angles[2] = acos( sum( cell_matrix[:,1].*cell_matrix[:,3] )/(length[1]*length[3]) )*tau
    angles[3] = acos( sum( cell_matrix[:,1].*cell_matrix[:,2] )/(length[1]*length[2]) )*tau
    #---------------------------------------------------------------------

    # Returns a cell parameter object
    return Cell_param( length, angles )
end
# Converts a cell trajectory in matrix form into a vector of cell params
# NB: We may want to create an additionnal object if not here, then in traj, to
# avoid dealing with an array of structure and have a structure of array instead
function matrix2Params( cell_matrices::Array{T1,3} ) where { T1 <: Real }
    # Arguments:
    # - cell_matrices: cell matrix in tensor form, time is assumed to be the third component of the matrix
    # Output:
    # - cell_params: Vector of Cell_param object that describes

    # Compute number of step
    nb_step = size( cell_matrices )[3]

    # Initialize output
    cells_params = Vector{ Cell_param }(undef, nb_step )

    # Loop over time
    for step=1:nb_step
        cells_params[step] = matrix2Params( cell_matrices[:,:,step] )
    end

    # Return object
    return cells_params
end
# Converts cell parameters into a cell matrix
function params2Matrix( cell_params::T1 ) where { T1 <: Cell_param }
    # Argument:
    # - cell_params: Cell_Param (see above), contains cell parameters (a,b,c),(alpha,beta,gamma)
    # Output:
    # - matrix: cell_matrix, contains all information about the cell in matrix form

    # Initialize output
    matrix=zeros(3,3)

    # Makes a copy of the parameters
    lengths=copy( cell_params.length)
    angles=copy(cell_params.angles)*pi/180 # Converts angles from degrees to radians

    # Compute matrix parameters
    #-----------------------------------------------------
    matrix[1,1] = lengths[1]
    matrix[1,2] = lengths[2]*cos( angles[3] )
    matrix[1,3] = lengths[3]*cos( angles[2] )
    matrix[2,2] = lengths[2]*sin( angles[3] )
    matrix[2,3] = lengths[3]*( cos( angles[1] ) - cos( angles[2] )*cos( angles[3] ) )/sin( angles[3] )
    volume=sqrt( 1 + 2*cos(angles[1])*cos(angles[2])*cos(angles[3]) -cos(angles[1])^2 -cos(angles[2])^2 -cos(angles[3])^2 )
    matrix[3,3] = lengths[3]*volume/sin( angles[3] )
    #-----------------------------------------------------

    # Return matrix
    return matrix
end
# Converts a trajectory of cell from a cell_params into a tensor form with time on the third dimension
function params2Matrix( cells_params::Vector{T1} ) where { T1 <: Cell_param }
    # Argument:
    # - cells_params: Vector of Cell_Param containing the trajectory of the cell in parameters form
    # Output:
    # - Tensor contaning the cell trajectory as cell matrices, with time on the third dimension of the tensor

    # Get number of step of trajectory
    nb_step = size(cells_params)[1]

    # Initialize output
    cells_ = zeros(Real, 3, 3, nb_step )

    # Loop over time
    #----------------------------
    for step=1:nb_step
        # Conversion
        cells_[:,:,step] = params2Matrix( cells_params[step] )
    end
    #----------------------------

    # Return tensor
    return cells_
end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Inverts the matrix cell, if possible
function invertCell( cell_matrix::Array{T1,2} ) where { T1 <: Real }
    # Argument
    # - cell_matrix: cell matrix describing the cell
    # Output
    # - inverse of the cell matrix (3x3 real scalar matrix) if inversion if possible, false otherwise

    # If det( matrix ) is not null then matrix is inversible (det comes from LinearAlgebra)
    if det( cell_matrix ) != 0
        # Computes and return inverse
        return LinearAlgebra.inv(cell_matrix)
    # Otherwise inversion is impossible, returns false
    else
        return false
    end
end
# Returns the inverse of the matrix cell corresponding to the cell parameters given as Cell_param
function invertCell( cell_params::T1 ) where { T1 <: Cell_param }
    # Argument
    # - cell_params: Cell_param describing the cell
    # Output
    # - 3x3 matrix containing the inverse of the cell matrix, or False if inversion is not possible

    # Computes the matrix linked to the params
    out = invertCell( params2Matrix(cell_params) )

    # If the inversion fails, return false
    if out == false
        return false
    # Otherwise returns the inverse of the cell matrix
    else
        return out
    end
end
#-------------------------------------------------------------------------------


# Computation of the volume of the cell
#-------------------------------------------------------------------------------
# Computes the volume using the cell matrix describing the cell
function getVolume( cell_matrix::Array{T1,2} ) where { T1 <: Real }
    # Argument:
    # - cell_matrix: cell matrix describing the target cell
    # Output:
    # - Volume of the cell (Real Scalar, strictly positive)

    # Computes and returns the volume
    return LinearAlgebra.det(cell_matrix)
end
# Computes the volume using the cell parameters
function getVolume( cell_param::T1 ) where { T1 <: Cell_param }
    # Argument
    # - cell_param: Cell_param object that contains the information about the cell
    # Output:
    # Volume of the cell (Real scalar, strictly positive)

    # Computes and return the volume of the cell
    return LinearAlgebra.det( params2Matrix( cell_param ) )
end
#-------------------------------------------------------------------------------

# Wrap functions (put the atoms back into the cell)
#-------------------------------------------------------------------------------
# Wrap a scalar position, within a box defined by a given length
function wrap( position::T1, length::T2 ) where { T1 <: Real, T2 <: Real}
    # Argument
    # - position: the position to wrap, ( real scalar, strictly positive )
    # - length:   the maximum allowed length for that position
    # Output:
    # - position (transformed): the position wrapped ( real scalar, strictly positive)

    # Wrap the position
    #-------------------------------------------
    # Determine
    sign=-1
    if position < 0
        sign=1
    end
    while position < 0 || position > length
        position = position + sign*length
    end
    #-------------------------------------------

    # Return the wrapped position
    return position
end
# Wrap atoms into a cell using a cell matrix for the cell
function wrap( atoms::T1, cell::Array{T2,2} ) where { T1 <: atom_mod.AtomList, T2 <: Real }
    # Arguments:
    # - atoms: atomic positions, in the format of AtomList
    # -  cell: cell matrix describing the cell to wrap the atom in
    # Output:
    # - atoms (transformed): atomic positions wrapped, in the format of AtomList

    # Getting Scaled positions
    atoms.positions = getTransformedPosition( atoms.positions, LinearAlgebra.inv( cell.matrix ) )

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

    # Returning the new positions as an AtomList
    return atoms
end
# Wrap atoms into a given cell using cell_param object for the cell
function wrap( atoms::T1, cell::T2 ) where { T1 <: atom_mod.AtomList, T2 <: Cell_param }
    # Arguments
    # - atoms: atomic positions in AtomList format
    # - cell : cell information
    # Output:
    # - atoms (transformed): atomic positions, wrapped in the cell, AtomList format

    # Loop over all the atoms contained in the frame
    for i=1:size(atoms.positions)[1]
        for j=1:3
            atoms.positions[i,j] = wrap( atoms.positions[i,j],cell.length[j])
        end
    end

    # Returning frame
    return atoms
end
# Wrap atoms in a set of atomic frames (AtomList format) into a given cell, using Cell_param vector for the cells
function wrap( traj::Vector{T1}, cell::T2 ) where { T1 <: atom_mod.AtomList, T2 <: Cell_param }
    # Arguments:
    # - traj : set of AtomList frames, includes atomic positions for each frame.
    # - cell : information regarding the cell
    # Output:
    # - traj (transformed): set of AtomList frames, positions wrapped in the cell

    # Get the number of frames
    nb_frames=size(traj)[1]

    # Loop over the frames
    for frame =1:nb_frame
        # Wrapping each frame
        traj[frame] = wrap( traj[frame], cell )
    end

    # Return the wrapped set
    return traj
end
# Wrap atoms in an AtomMolList frame into a given cell, using Cell_param vector for the cell
function wrap( molecules::T1, cell::T2 ) where { T1 <: atom_mod.AtomMolList, T2 <: Cell_param }
    # Arguments:
    # - molecules: AtomMolList frame containing atomic positions
    # - cell: Cell_param describing the cell
    # Output:
    # - molecules (transformed): Same frame but with wrapped atomic positions

    # Loop over frames
    for atom=1:size(molecules.positions)[1]
        for j=1:3
            molecules.positions[i,j] = wrap( molecules.positions[atom,j], cell.length[j] )
        end
    end

    # Return the same set, but with position wrapped
    return molecules
end
# Wrap atomic positions in a given cell, using an array for the positions the Cell_param structure for the cell
function wrap( positions::Array{T1,2}, cell::T2 ) where { T1 <: Real, T2 <: Cell_param }
    # Arguments
    # - positions : array containing atomic positions
    # - cell      : Cell_param describing the cell
    # Output
    # - positions (transformed): array containing the atomic positions (wrapped)

    # Loop over atoms
    for atom=1:size(positions)[1]
        # Loop over dimension
        for i=1:3
            # Wrapping in each dimension
            positions[atom,i] = wrap( positions[atom,i],cell.length[i] )
        end
    end

    # Returning the wrapped positions
    return positions
end
# Wrapping atomic position of a single atom (as a vector), in a cell gfiven by a Cell_param Object
function wrap( positions::Vector{T1}, cell::T2 ) where { T1 <: Real, T2 <: Cell_param }
    # Argument:
    # - positions : the atomic position of the atom (vector of real numbers)
    # - cell : cell described by a Cell_param

    # Loop over dimensions
    for i=1:3
        # Wrapping over dimensions
        positions[i] = wrap( positions[i],cell.length[i] )
    end

    # Return the atomic position
    return positions
end
#-------------------------------------------------------------------------------

# Unwrapping functions for molecule analysis
#---------------------------------------------------------------------------------------------
# Unwrap a target atom using a reference atom, all positions are in array, cell given as a Cell_param
# - Works for Orthorombic cell
# - Modifies Array
function unWrapStructureOrtho!( positions::Array{T1,2}, reference::T2, target::T3, cell::T4  ) where { T1 <: Real, T2 <: Int, T3 <: Int, T4 <: Cell_param }
    # Arguments:
    # - positions: trajectory of a single atom in matrix form (n_step,3) (real)
    # - reference: reference atom (int)
    # - target: target atom to unwrap (int)
    # - cell: Cell_param that describes the cell
    # Output:
    # - Nothing, modifications made direclty on array

    # Loop over the dimensions
    for i=1:3
        # Compute distance between target and reference
        dist = ( positions[reference,i] - positions[target,i] )
        # If distance is larger than half the size of the box, unwrap
        if dist > cell.length[i]*0.5
            positions[target,i] += cell.length[i]
        # If distance is smaller than - half the size of the box, unwrap
        elseif dist < -cell.length[i]*0.5
            positions[target,i] -= cell.length[i]
        end
    end

    return
end
# Unwrap a target atom using a reference atom, positions in AtomList, cell given in Cell_param
# - Works for orthorombic cell only
# - Modifies AtomList
function unWrapStructureOrtho!( structure::T1, reference::T2, target::T3, cell::T4  ) where { T1 <: atom_mod.AtomList, T2 <: Int, T3 <: Int, T4 <: Cell_param }
    # Argument
    # - structure: AtomList that contains the atomic positions
    # - origin: index of the reference atom
    # - target: index of the target atom
    # - cell: Cell_param that describe the cell
    # Output
    # - Nothing, modifications are made directly on structure

    # Loop over dimensions
    for i=1:3
        # Compute distance between reference and target
        dist = ( structure.positions[reference,i] - structure.positions[target,i] )

        # If distance is larger than half the cell, unwrap
        if dist > cell.length[i]*0.5
            structure.positions[target,i] += cell.length[i]
        # If distance is smaller than minus half the cell, unwrap
        elseif dist < -cell.length[i]*0.5
            structure.positions[target,i] -= cell.length[i]
        end
    end

    return
end
# Unwrap all atoms within a single molecule, to generate the actual form that is broken by PBC
# Recursive (may generate infinite loop)
# Works for orthorombic cell only
function unWrapStructureOrthoOnce!( visited::Vector{T1}, adjacency_table::Vector{T2}, positions::Array{T3,2} , cell::T4, target::T5, index_atoms::Vector{T6} ) where { T1 <: Int, T2 <: Any, T3 <: Real, T4 <: Cell_param, T5 <: Int, T6 <: Int }
    # Arguments
    # - visited: map of all atoms that were previously visited by the function
    # - adjacency_table: matrix describing which atoms are bonded (int: 0=non bonded, 1=bonded)
    # - positions: atomic positions in array form
    # - cell: Cell_param describing cell
    # - index: conversion table that matchesthe neighbor number to actual index atom (int, vector)
    # Output:
    # Nothing, modifications active on positions

    # Target atom is marked as visited
    visited[target]=1

    # Get number of neighbor of target
    nb_neighbor=size(adjacency_table[target])[1]

    # Loop over neighbors
    for neigh=1:nb_neighbor
        # If neighbor was already visited, move to the next one
        if visited[ adjacency_table[target][neigh] ] == 0
            # Unwrap the neighbor that was not visited
            unWrapStructureOrtho!( positions, index_atoms[target], index_atoms[ adjacency_table[target][neigh] ], cell )
            # Goes to unwrap all the neighbors of the neighbor atom
            unWrapStructureOrthoOnce!( visited, adjacency_table, positions, cell, adjacency_table[target][neigh], index_atoms )
            # When all neighbors of the branch were visited and unwrap, we go to the next one.
        end
    end

    return
end
#-------------------------------------------------------------------------------

# Functions dealing with Cartesian <-> Reduced coordinates
# NB: TEST THOSE FUNCTIONS
#-------------------------------------------------------------------------------
# Transforms the position of an atom (real vector) using the matrix given in parameter
# - From cartesian to reduced with the inverse cell matrix
# - From reduced to cartesian coordinates with the cell matrix
function getTransformedPosition( target_vector::Vector{T1}, cell_matrix::Array{T2,2} ) where {T1 <: Real, T2 <: Real }
    # Argument:
    # - target_vector: original position (real vector of size 3)
    # - cell_matrix: 3x3 matrix, either the cell matrix or its inverse
    # Output:
    # - vector: transformed position (real vector of size 3)

    # Initialize output
    vector=zeros(3)

    # Loop over dimension 1
    for i=1:3
        # Loop over dimension 2
        for j=1:3
            vector[i] += cell_matrix[i,j]*target_vector[j]
        end
    end

    # Returns the transformed atomic position
    return vector
end
# Transforms positions of atoms (real matrix form) using the matrix given in parameter
# - From cartesian to reduced with the inverse cell matrix
# - From reduced to cartesian coordinates with the cell matrix
function getTransformedPosition( target_matrix::Array{T1,2}, cell_matrix::Array{T2,2} ) where {T1 <: Real, T2 <: Real }
    # Argument
    # - target_matrix: atomic position in real matrix form (nb_atoms,3)
    # - cell_matrix: either cell matrix or its inverse, (3x3 real matrix)
    # Output
    # - matrix_transformed: Transformed positions in matrix form (nb_atoms,3)

    # Initialize output
    matrix_transformed = zeros( nb_point, 3 )

    # Loop over atoms
    for atom=1:size(target_matrix)[1]
        # Transforms each atomic positions
        matrix_transformed[ atom, : ] = getTransformedPosition( target_matrix[ atom,: ], cell_matrix )
    end

    # Returns the transformed positions
    return matrix_transformed
end
# Transforms cartesian position of a single atom into reduced coordinates for an atom using matrix
function cartesian2Reduced( vector::Vector{T1}, cell_matrix::Array{T2,2} ) where { T1 <: Real , T2 <: Real }
    # Arguments
    # - vector: atomic position cartesian coordinates (real vector, size=3)
    # - cell_matrix: cell matrix describing the cell
    # Output:
    # - vector containing reduced coordinates of the atom (real vector, size=)

    # Compute the inverse of the cell matrix
    cell_mat_inverse = invertCell( cell_matrix )

    # If cell matrix is not inversible, returns false
    if cell_mat_inverse == false
        return false
    end

    # Compute and returns the reduced coordinates of the atom
    return getTransformedPosition( vector, cell_mat_inverse )
end
# Transforms cartesian position of a set of atoms (matrix form) into reduced coordinates for an atom using matrix
function cartesian2Reduced( target_matrix::Array{T1,2}, cell_matrix::Array{T2,2} ) where { T1 <: Real, T2 <: Real }
    # Arguments:
    # - target_matrix: Contains the atomic positions as a real matrix (nb_atoms,3) in cartesian coordinates
    # - cell_matrix: 3x3 matrix describing the cell
    # Output:
    # - scaled_positions: transformed positions into reduced coordinates (real matrix (nb_atoms,3))
    # if matrix is not inversible, returns false

    # Computes the inverse of the cell matrix
    inv_cell_matrix=invertCell(cell_matrix)

    # If cell_matrix is not inversible, returns false
    if inv_cell_matrix == false
        return false
    end

    # Initialize output (may not be necessary)
    reduced_positions=copy(target_matrix)

    # Loop over all atoms
    for atom=1:size(target_matrix)[1]
        # Transform each atomic position into reduced form
        reduced_positions[atom,:] = getTransformedPosition( target_matrix[atom,:], inv_cell_matrix )
    end

    # Return reduced coordinates of all atoms
    return reduced_positions
end
# Transforms reduced coordinates of an atom into cartesian coordinates using cell matrix for the cell
function reduced2Cartesian( positions_reduced::Vector{T1}, cell_matrix::Array{T2,2} ) where { T1 <: Real, T2 <: Real }
    # Argument
    # - positions_reduced: reduced coordinates of a single atom
    # - cell_matrix: matrix that describes the cell
    # Output:
    # - new_positions: cartesian coordinates of the atom corresponding to the reduced

    # Initialize output
    new_positions=zeros(Real,3)

    # Loop over cartesian coordinates dimension
    for i=1:3
        # Loop over reduced coordinates dimension
        for j=1:3
            new_positions[i] += positions_reduced[j]*cell_matrix[i,j]
        end
    end

    # Returns the cartesian coordinates
    return new_positions
end
# Transforms reduced coordinates of a set of atoms (real matrix form) into cartesian coordinates, using the cell matrix
function reduced2Cartesian( positions_reduced::Array{T1,2}, cell_matrix::Array{T2,2} ) where { T1 <: Real, T2 <: Real }
    # Arguments:
    # - positions_reduced: atomic position in reduced coordinates, matrix form (nb_atoms,3)
    # - cell_matrix: matrix describing the cell
    # Output:
    # - positions_reduced, transformed: the atomic positions in cartesian coordinates

    # Loop over all the atoms
    for atom=1:size(positions_reduced)[1]
        # Transforms the positions to cartesian
        positions_reduced[atom,:] = reduced2Cartesian( positions_reduced[atom,:], cell_matrix )
    end

    # Returns the cartesian coordinates
    return positions_reduced
end
#-------------------------------------------------------------------------------

# Computation of distances between atoms
# TODO Check and sort this
#-------------------------------------------------------------------------------
# Compute distance in 1 dimension
function distsq1D( x1::T1, x2::T2, a::T3 ) where { T1 <: Real, T2 <: Real, T3 <: Real }
    # Arguments:
    # - x1, x2: atomic positions
    # - a: length of the cell
    # Output:
    # squared distance between x1 and x2 with PBC

    # Compute distance
    dx=x1-x2

    # Account for PBC
    if dx > a*0.5
        dx -= a
    end
    if dx < -a*0.5
        dx += a
    end

    # Return squared distance
    return dx*dx
end
# Computes the distance between two atoms
# - works for orthorombic only
function distanceOrtho( v1::Vector{T1}, v2::Vector{T2}, cell_length::Vector{T3} ) where { T1 <: Real, T2 <: Real, T3 <: Real }
    # Argument
    # - v1, v2: atomic positions, real vectors of size 3
    # - cell_length: lengths of the cell in all three dimension
    # Output
    # - distance between v1 and v2 (real value, positive)

    # Initialize output
    dist=0

    # Loop over dimensions
    for i=1:size(v1)[1]
        # Computes squared distance in each dimension
        dist += distsq1D( v1[i], v2[i], cell_length[i] )
    end

    # Return distance
    return sqrt(dist)
end
# Computes the distance between two atoms in a set of positions (matrix form)
# - works for orthorombic only
function distanceOrtho( positions::Array{T1,2}, cell::T2, index1::T3, index2::T4 ) where { T1 <: Real,  T2 <: Cell_param, T3 <: Real, T4 <: Real }
    # Argument
    # - positions: matrix containing positions of all atoms (nb_atoms,3)
    # - cell: Cell_param containing all info about the cell
    # - index1, index2: index of target atoms (int)
    # Output:
    # - dist: distance between atom1 and atom2 (real positive scalar)

    # Initialize output
    dist=0

    # Loop over dimensions
    for i=1:3
        # Compute distance in each dimension
        dist +=  dist1D( positions[index1,i], positions[index2,i], cell.length[i] )
    end

    # Returns distance
    return sqrt(dist)
end
# Computes the distance between two atoms in a set of positions (AtomList form)
# - works for orthorombic only
function distance( atoms::T1, cell::T2, index1::T3, index2::T4 ) where { T1 <: atom_mod.AtomList, T2 <: Cell_param , T3 <: Int, T4 <: Int }
    # Arguments
    # - atoms: AtomList that contains all atomic positions
    # - cell: Cell_param that contains all info about the cell
    # - index1, index2: indexes of the two target atom (int)
    # Output
    # - dis: distance between two atoms index1 and index2

    # Initialize output
    dis=0

    # Loop over dimensions
    for i=1:3
        # Compute distance squared in each dimensions
        dis += dist1D( atoms.positions[index1,i],atoms.positions[index2,i], cell.length[i] )
    end

    # Returns the distance between two atoms
    return sqrt(dis)
end
# Computes the distance between two atoms in a set of positions (AtomList form)
# - works for orthorombic only
# - can Wrap atoms before computing distance
function distanceOrtho( atoms::T1, cell::T2, index1::T3, index2::T3, wrap::T4 ) where { T1 <: atom_mod.AtomList, T2 <: Cell_param,  T3 <: Int, T4 <: Bool }
    # Arguments
    # - atoms: AtomList that contains atomic positions between atoms
    # - cell: Cell_param that contains all information about the cell
    # - index1, index2: indexes of the target atoms (int)
    # - wrap: whether or not to wrap atoms
    # Output
    # a real vector, the distance between atoms index1 and index2

    # Wrap (or not atoms)
    if (  wrap )
        wrap( atoms, cell )
    end

    # Computes and returns the distance between atoms index1 and index2
    return distance( atoms, cell, index1, index2 )
end
# Computes the distance between two atoms in a set of positions (AtomMolList form)
# - works for orthorombic only
function distanceOrtho( atoms::T1, cell::T2, index1::T3, index2::T4 ) where { T1 <: atom_mod.AtomMolList, T2 <: Cell_param , T3 <: Int , T4 <: Int }
    # Argument
    # - atoms: AtomMolList containing all atomic positions
    # - cell: Cell_param describing the cell
    # - index1, index2 : indexes of the target atoms (int)
    # Output
    # dist: contains the distance between atoms index1 and index2

    # Initialize output
    dist=0

    # Loop over dimension
    for i=1:3
        # Compute distance squared in all dimensions
        dist += dist1D( atoms.positions[index1,i], atoms.positions[index2,i], cell.length[i] )
    end

    # Return distance
    return sqrt(dist)
end
# Computes distance in Reduced Coordinates
function distanceReduced( v1_scaled::Vector{T1}, v2_scaled::Vector{T2}, cell_matrix::Array{T3,2} ) where { T1 <: Real, T2 <: Real, T3 <: Real }
    # Arguments
    # - v1_scaled, v2_scaled: position in reduced coordinates
    # - cell_matrix: cell matrix describing the cell
    # Output:
    # - distance in reduced coordinates between v1_scaled and v2_scaled

    # Initialize output
    ds=zeros(3)

    # Distance + Min Image Convention
    for i=1:3
        ds[i] = v1_scaled[i] - v2_scaled[i]
        ds[i] = ds[i] - round(ds[i])
    end

    # Descaling distance vector
    ds = getTransformedPosition( ds, cell_matrix )

    # Returns the distance in reduced space
    return sqrt( dot( ds, ds ) )
end
# Computes the distance between two atoms, works for all cells (but slower)
function distance( v1::Vector{T1}, v2::Vector{T2}, cell_matrix::Array{T3,2} ) where { T1 <: Real, T2 <: Real, T3 <: Real }
    # Argument
    # - v1, v2: positions in reduced space
    # - cell_matrix: cell matrix describing the cell
    # Output
    # - distance between v1 and v2, real positive scalar

    # Converts positions from cartesian to reduced
    v1_scaled=getScaledPosition(v1,cell_matrix)
    v2_scaled=getScaledPosition(v2,cell_matrix)

    # Initialize distance in scaled space
    ds=zeros(3)

    # Compute Distance in scaled space with PBC
    for i=1:3
        ds[i] = v1_scaled[i]-v2_scaled[i]
        ds[i] = ds[i] - round(ds[i])
    end

    # Descaling distance
    ds=getTransformedPosition(ds,cell_matrix)

    # Returns distance
    return sqrt(dot(ds,ds))
end
# Computes the distance between two atoms, works for all cells (but slower) - redundant
function distance( v1::Vector{T1}, v2::Vector{T2}, cell_matrix::Array{T3,2}, inv_cell_matrix::Array{T4,2} ) where { T1 <: Real, T2 <: Real, T3 <: Real, T4 <: Real }

    # Converts cartesian coordinates to reduced coordinates
    v1_scaled=getTransformedPosition(v1,inv_cell_matrix)
    v2_scaled=getTransformedPosition(v2,inv_cell_matrix)

    # Initialize distance in reduced coordinates
    ds=zeros(3)

    # Computes Distance with PBC
    for i=1:3
        ds[i] = v1_scaled[i]-v2_scaled[i]
        ds[i] = ds[i] - round(ds[i])
    end

    # Descaling of distance vector
    ds=getTransformedPosition(ds,cell_matrix)

    # Returns distance
    return sqrt(dot(ds,ds))
end
# Computes the distance between two atoms, works for all cells (but slower), with Cell_param to describe the cell
function distance( v1::Vector{T1}, v2::Vector{T2}, cell_params::T3 ) where { T1 <: Real, T2 <: Real, T3 <: Cell_param }
    # Arguments
    # - v1, v2: atomic positions of two atoms
    # - cell_params: Cell_param object that describes the cell
    # Output
    # - Distance between v1 and v2

    # Computes the matrix related to the parameters,
    # Computes the distance
    # Returns the distance
    return distance( v1, v2, params2Matrix( cell_params ) )
end
#---------------------------------------------------------------------------

# Computes the velocities in a traj by finite element method from the positions
# NB: This is not a perfect method, use only with small timestep (a few at most)
#---------------------------------------------------------------------------
function velocityFromPosition( traj::Vector{T1}, dt::T2, dx=1::T3 ) where { T1 <: atom_mod.AtomList, T2 <: Real, T3 <: Real }
    # Arguments:
    # traj: trajectory in the form of a vector of AtomList, contains all atomic positions
    # dt: timestep for the trajectory
    # dx: conversion of the distance units (optional)
    # Output:
    # velocities: array(nb_step,nb_atoms,3) containing velocities of atoms

    # Initialize output
    velocities = zeros( size(traj)[1]-1, size(traj[1].names)[1], 3 )

    # Loop over all step (minus 1, as we can't compute the velocity for the last step)
    for step=1:size(traj)[1]-1
        # Loop over the atoms
        for atom=1:size(traj[1].names)[1]
            # Loop over dimensions
            for i=1:3
                # Compute velocity component i for atom at step
                velocities[step,atom,i] = dx*( traj[step].positions[atom,i] - traj[step+1].positions[atom,i] )/dt
            end
        end
    end

    # Return velocities
    return velocities
end
#---------------------------------------------------------------------------

# Compressing cell function
# NB: Not sure it works, and will create serious issues for molecular systems
# in already reduced states
#---------------------------------------------------------------------------
# Compress the cell using cell_param by a given factor
# - works for orthorombic cell
function compressParams( cell::T1, fracs::Vector{T2} ) where { T1 <: Cell_param, T2 <: Real }
    # Arguments
    # - cell: Cell_param describing the original cell
    # - fracs: vector of size 3 containing the factor by which to multiply the lengths of the cell (real,>0)
    # Output
    # - cell: cell_param containing the information of the new cell

    # Compression of the cell by reducing its lengths
    for i=1:3
        cell.lengths[i] *= fracs[i]
    end

    # Returns the new, compressed, cell
    return cell
end
#---------------------------------------------------------------------------

# Checks for infinite molecules
#---------------------------------------------------------------------------
# Context: used to determine whether molecules where looping over themselves
# through the PBC in some cases where they form through polymerization.
#TODO: Check that they actually work
#TODO: Check which of these is actually useful and which isn't
#---------------------------------------------------------------------------
# Determine whether the molecule to which the target atom belong is infinite (self-looping) - graph exploration method
# - Recursive, careful with it
# NB: This should be tested
function isInfiniteChain( visited::Vector{T1}, matrix::Array{T2,2}, adjacency_table::Vector{T3}, positions::Array{T4,2} , cell::T5, target::T6, index_atoms::Vector{T7}, cut_off::T8 ) where { T1 <: Int, T2 <: Real, T3 <: Any, T4 <: Real, T5 <: cell_mod.Cell_param, T6 <: Int, T7 <: Int, T8 <: Real }
    # Arguments
    # - visited: vector (int) containing the information about whether or not the atom was already visited by the exploration (1 visited, 0 no visited)
    # - matrix: distance matrix between all atoms in the set  (real: (nb_atoms,nb_atoms))
    # - adjacency_table: connection matrix, determining whether atoms are bonded or not (1=bonded,0=not bonded)
    # - positions: matrix containing atomic positions (nb_atoms,3) (real)
    # - cell: Cell_param describing the cell
    # - target: target atom to start with (int)
    # - index_atoms: bookeeps the actual index of all atoms (vector of int)
    # - cut_off: cut_off to determine bonding (real)
    # Output
    # - Is the molecule an infinite chain? (Bool)
    # - Was the chain exploration went ok? (Bool)

    # How this works:
    # This function tries to unwrap progressively a tree-graph corresponding to the molecule
    # to which the target atom belongs to. If by unwraping progressively (and only once per atom), a bond disapear,
    # this means the molecule loops over itself and is infinite

    # Checking the target atom as visited
    visited[target] = 1

    # Getting the number of neighbor of target atom
    nb_neighbor = size( adjacency_table[target] )[1]

    # Loop over all neighbors
    for neigh=1:nb_neighbor

        # Checking that neighbor is not already visited
        if visited[adjacency_table[target][neigh]] == 0

            # Unwrapping neighbor of target
            unWrapOrtho!( positions, index_atoms[target], index_atoms[ adjacency_table[target][neigh] ], cell )

            # Checking whether atoms are appart (when they should not be)
            if geom.distance( positions[ index_atoms[target], : ], positions[ index_atoms[ adjacency_table[target][neigh] ] , : ] ) > cut_off
                return true, true
            end

            # Sending the search to the neighbor
            isinf, isok = isInfiniteChain(visited,matrix,adjacency_table,positions,cell,adjacency_table[target][neigh],index_atoms,cut_off)

            # Something went wrong
            if ! isok
                return true, false
            end

            # Chain is infinite
            if ! isinf
                return true, true
            end

        # If the neighbor was visited but that the distance is too large with regard to the cut-off, infinite chain spotted
        elseif geom.distance(positions[index_atoms[target],:],positions[ index_atoms[adjacency_table[target][neigh]] ,: ] ) > cut_off
            # Spotted infinite loop; stops the search
            return true, true
        end
    end

    # If we get here, the molecule was properly reconstructed and is not infinite
    return false, true
end
# Check if molecule is infinite by unwraping once all atoms following a graph exploration method
function checkInfiniteChain( matrix::Array{T1,2},  positions::Array{T2,2} , cell::T3, molecule_indexs::Vector{T4}, cut_off::T5 ) where { T1 <: Real, T2 <: Real, T3 <: Cell_param, T4 <: Int, T5 <: Real }
    # Arguments:
    # - matrix: adjacency matrix, 1=atoms i,j are bonded, 0 = atoms i,j are not bonded
    # - positions: Array contaning the positions of the atoms (nb_atoms,3)
    # - cell :  Cell_param containing all info on the cell
    # - molecule_indexs: index of all atoms in the molecules
    # - cut_off : cut-off used to determine bonding in the adjacency matrix
    # Output
    # - Is the molecule infinite? (Bool)

    # Compute neighbors of all atoms
    adjacent_molecule = graph.getAllAdjacentVertex( matrix )

    # Compute the size of the molecule
    size_molecule = size( matrix )[1]

    # Creat a vector to bookkeep which atoms were visited
    visited=zeros(Int,size_molecule)

    # Unwrapping all atoms but only once
    cell_mod.unWrapOrthoOnce( visited, matrix , adjacent_molecule, positions, cell, 1, molecule_indexs )

    # Loop over atoms of the molecule
    for atom=1:size_molecule
        # Loop over neighbors
        for adj=1:size(adjacent_molecule[atom])[1]
            # Checking bonds correspondance
            if norm( positions[ molecule_indexs[atom], :] - positions[ molecule_indexs[adjacent_molecule[atom][adj]], :] ) > cut_off
                # If some atoms are no longer bonded, then we have an infinite chain
                return true
            end
        end
    end

    # If we get there, the molecule isn't infinite
    return false
end
# Find the atoms that are bonded but whose bonds are cut by the one unwrap policty within a single molecule
function findUnlinked( matrix::Array{T1,2},  positions::Array{T2,2} , cell::T3, molecule_indexs::Vector{T4}, cut_off::T5  ) where { T1 <: Real, T2 <: Real, T3 <: Cell_param, T4 <: Int, T5 <: Real }
    # Arguments
    # - matrix: adjacency matrix of the molecule, if (i,j) = 1 atoms i,j are bonded, if 0 they are not
    # - positions: real array (nb_atoms,3) contains the positions of atoms
    # - cell: Cell_param describing the cell
    # - molecule_indexs: bookeep the indexes of the atoms of the molecule
    # - cut_off: cut-off used to create the adjacency matrix of the molecule (real scalar > 0 )
    # Output
    # - list: Array (nb_pair,2) that contains all pairs of atoms that are unlinked by unwrapping

    # Init list of unlinked atoms
    list=Array{Int}(undef,0,2)

    # Get the size of the molecule
    size_molecule=size(matrix)[1]

    # Loop over atoms in the molecule
    for atom=1:size_molecule
        # Second loop over atoms of the molecule
        for atom2=atom+1:size_molecule
            # Checking that atoms that should be bonded are bonded
            # matrix[atom1,atom2] : the truth value
            # norm(...) : whether atoms are bonded post unwrapping
            if matrix[atom,atom2] == 1 && norm(positions[atom,:]-positions[atom2,:]) > cut_off
                # If there is no unbonded atoms
                if size(list)[1] == 0
                    # Adding atoms that are problematic to the list
                    list = vcat( list, [ molecule_indexs[atom] molecule_indexs[atom2] ] )
                else
                    # Checking that we haven't already added the pair
                    check = true
                    # Loop over all atoms in the list
                    for i=1:size(list)[1]
                        if list[i,:] == [ molecule_indexs[atom], molecule_indexs[atom2] ]
                            check = false
                        end
                    end
                    # If the atoms are not in the list...
                    if check
                        # We add the atom to the list
                        list = vcat( list, [ molecule_indexs[atom] molecule_indexs[atom2] ] )
                    end
                end
            end
        end
    end

    # Returns the list of unlinked atoms if any
    return list
end
#---------------------------------------------------------------------------

# Supercell building functions
#-------------------------------------------------------------------------------
# Computes how atoms of a unit cell must move when growing for a supercell
function computeMoveVector( index::Vector{T1}, cell_matrix::Array{T2,2} ) where { T1 <: Int, T2 <: Real }
    # Arguments
    # - index: direction of the growth of the original size (vector, size 3, int)
    # - cell_matrix: cell matrix, 3x3 real matrix describing the cell
    # Output
    # - moveVector: A vector by which to move all atoms of a cell to create a supercell

    # Initialize output
    moveVector = zeros(Real,3)

    # Loop over dimension
    for i=1:3
        # Second loop over dimension
        for j=1:3
            # Compute the move
            moveVector[i] += index[i]*cell_matrix[i,j]
        end
    end

    # Return the move vector
    return moveVector
end
# Growing the cell into a supercell using n_grow indexes, cell described by cell matrix
function growCell( cell::Array{T1,2}, n_grow::Vector{T2} ) where { T1 <: Real, T2 <: Int }
    # Argument
    # - cell : cell matrix of the original cell
    # - n_grow: index by which to grow the cell in all direction (vector, 3, int)
    # Output
    # - cell2: cell matrix of the new cell

    # Initialize the output cell by copying the original
    cell2 = copy(cell)

    # Constructing new cell
    for i=1:3
        cell2[:,i] = cell[:,i]*n_grow[i]
    end

    # Returns the modified cell
    return Cell_matrix( cell2 )
end
# Growing the cell into a supercell using n_grow indexes, cell described by cell_param
function growCell( cell::T1, n_grow::Vector{T2} ) where { T1 <: Cell_param, T2 <: Int }
    # Arguments
    # - cell: Cell_param of the original cell
    # - n_grow: indexes used to grow the cell in all 3 dimensions
    # Output
    # - cell matrix modified for the supercell

    # Computes the cell matrix,
    # Compute and returns the supercell
    return growCell( params2Matrix(cell), n_grow )
end
# Making a supercell using the n_grow indexes to decide in which direction to grow
function makeSuperCell( atoms::T1, cell_matrix::Array{T2,2}, n_grow::Vector{T3} ) where { T1 <: AtomList, T2 <: Real, T3 <: Int }
    # Get number of atoms in the original cell
    nb_atoms_base = size(atoms.names)[1]

    # Compute the number off atoms in the final cell
    nb_atoms_new = nb_atoms_base
    for i=1:3
        nb_atoms_new *= n_grow[i]
    end

    # Initialize final AtomList size
    new_atoms = AtomList( nb_atoms_new )

    # Loop over x direction
    for dirx=0:n_grow[1]-1
        # Loop over y direction
        for diry=0:n_grow[2]-1
            # Loop over z direction
            for dirz=0:n_grow[3]-1

                # Compute the direction in the x,y,z direction
                move_box = [dirx,diry,dirz]

                # Compute the move vector for all atoms of the original cell
                moveVector = zeros(Real,3)
                for xyz=1:3
                    for v123=1:3
                        moveVector[xyz] += move_box[v123]*cell_matrix[xyz,v123]
                    end
                end

                # Move all atoms by the move vector to populate the larger cell
                for atom = 1:nb_atoms_base
                    new_atoms.positions[count_,:] = atoms.positions[atom,:] .+ moveVector[:]
                    new_atoms.names[count_] = atoms.names[atom]
                    new_atoms.index[count_] = 1
                end

            end
        end
    end

    # Creating the supercell
    super_cell = growCell( cell_matrix, n_grow )

    # Sort Atoms by Z
    atom_mod.sortAtomsByZ!(new_atoms)

    # Changes the indexes to default
    for i=1:nb_atoms_new
        new_atoms.index[i] = i
    end

    return new_atoms, cell_mod.cellMatrix2Params( super_cell )
end
#-------------------------------------------------------------------------------

# Transforming cell from unorthorombic to orthorombic using a cut
# NB: Does not work for all structures, be extremely cautious with the results
#-------------------------------------------------------------------------------
function toOrthoByCut( cell_matrix::Array{T1,2} ) where { T1 <: Real }
    # Argument
    # - cell_matrix : original cell matrix
    # Output
    # - Cell_param that describes the newly created orthorombic cell

    # Initialize lengths of cell
    lengths = zeros( Real, 3 )

    # Aligne direction 1 along the x axis
    lengths[1] = LinearAlgebra.norm( cell_matrix[:,1] )

    # Loop to compute the lengths in the other directions
    for i=2:3
        for j=1:3
            lengths[j] += cell.matrix[i,j]
        end
    end

    # Returns the new cell (Cell_param)
    return cell_mod.Cell_param( lengths )
end
#-------------------------------------------------------------------------------

end
