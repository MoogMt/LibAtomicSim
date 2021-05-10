module cube_mod

# Loading necessary modules from LibAtomSim
using atom_mod
using cell_mod
using geom
using conversion
using periodicTable

# Loading necessary modules from Julia default repository
using LinearAlgebra

# Exporting functions
export readCube, writeCube
export applyPBCCube
export moveValues
export distanceSimpleCube
export carveCube
export getClosestIndex, computeDisplacementOrigin
export dataInTheMiddleWME, traceLine

# Volume structure to manipulate info
#-----------------------------------------------------------------------
mutable struct Volume
    # Variables
    #--------------------------------
    matrix::Array{Real}   # Tensor with cube data
    vox_vec::Array{Real}  # Vector describing the orientation of voxel vectors
    nb_vox::Vector{Int}   # Number of voxel per vector direction
    origin::Vector{Real}  # Shift between origin and location
    #--------------------------------

    # Constructors
    #---------------------------------------------------------------------------
    # Default constructor
    function Volume()
        new( Array{ Real, 3 }(undef, 1, 1, 1 ), Array{Real,2}(undef, 3, 3 ), Array{Int,1}(undef, 3 ),Array{Real,1}(undef, 3 ) )
    end
    # Cubic volume constructor
    function Volume( nb_vox::T1 ) where {T1 <: Real}
        # Argument
        # - nb_vox: number of voxels in all three direction (cubic cell)

        # Check that the number is positive
        if nb_vox > 0
            new( Array{Real}( nb_vox, nb_vox, nb_vox ), Array{Real}(3,3), [ nb_vox, nb_vox, nb_vox ] )
        # If not, sends an error message and return false
        else
            print("/!\\ Error: Wrong size for the number of voxels.\n")
            return false
        end
    end
    # Default constructor with all information specified
    function Volume( matrix::Array{T1,3}, vector_vox::Array{T2,2}, nb_vox::Vector{T3}, origin::Vector{T4} )  where { T1 <: Real, T2 <: Real, T3 <: Int, T4 <: Real }
        # Argument
        # - matrix: tensor (Real, nb_vox_x, nb_vox_y, nb_vox_z) containing all data
        # - vector_vox: voxel vectors describing the structure of the voxels (basically the reduced cell matrix, reduced by number of voxels in each vector direction)
        # - nb_vox: number of voxels in each dimensions
        # - origin: shift of the voxel grid compared to the position origin

        # Creates object
        new( matrix, vector_vox, nb_vox, origin )
    end
    #---------------------------------------------------------------------------
end
#-----------------------------------------------------------------------

# Reads a cube file and returns all information in the form of a volume
#-----------------------------------------------------------------------------------
function readCube( file_path::T1 ) where { T1 <: AbstractString }
    # Argument
    # - file_path: path to the input file
    # Output
    # - atom_list: AtomList object with atom positions, names and index
    # - cell_matrix: cell matrix describing cell
    # - volume: volume with ELF, density, or whatever is contained in the cube file

    # Max number of column for data
    nb_col = 6

    # Open input file
    handle_in = open( file_path )

    # Reads entire files and put all lines in tuple
    lines = readlines( handle_in )

    # Closing input file
    close( handle_in )

    # Ignoring first 2 lines
    offset=3

    # Number of atoms
    nb_atoms = parse(Int, split( lines[ offset ] )[1] );

    # Origin position of the density
    center=zeros(3)
    for i=1:3
        center[i] = parse(Float64, split( lines[ offset ] )[ i + 1 ])*conversion.bohr2Ang
    end

    # Number of voxels in each direction
    nb_vox=zeros(Int,3)
    for i=1:3
        nb_vox[i] = parse(Float64, split( lines[ i + offset ] )[1] )
    end

    # Parse Cell information
    cell_matrix=zeros(Real, 3, 3 )
    for i=1:3
        for j=1:3
            cell_matrix[ i, j ] = parse(Float64, split( lines[ i + offset ] )[ j + 1 ] )*conversion.bohr2Ang
        end
    end

    # Increase offset
    offset = 6

    # Parse Atomic Information
    atom_list=atom_mod.AtomList(nb_atoms);

    # Loop over atoms
    for atom=1:nb_atoms
        # Get atom names
        atom_list.names[atom] = periodicTable.z2Names( parse(Int, split( lines[ atom + offset ] )[1] ) )

        # Loop over dimensions
        for j=1:3
            # Parsing
            atom_list.positions[ j, atom ] = parse(Float64, split( lines[ atom + offset ])[ j + 2 ] )*conversion.bohr2Ang
        end
    end

    # Parse Volume Data
    nb_tot = nb_vox[1]*nb_vox[2]*nb_vox[3]

    # Init matrix for data
    matrix = zeros(Real, nb_vox[1], nb_vox[2], nb_vox[3] )

    # Update offset
    offset = round(Int, offset + nb_atoms + 1 )

    # Init position of matrix
    x=1; y=1; z=1;

    # Loop over volume lines
    for i=0:nb_tot/nb_col-1
        # Loop over columns
        for j=1:nb_col
            # index of current line
            index_ = round(Int, offset + i  )
            # volume data
            keys = split( lines[ index_ ] )
            #
            if size( keys )[1] < j
                break
            end
            # Storing data in matrix
            matrix[ x, y, z ] = parse(Float64, keys[ j ] )
            # Increment z
            z = z + 1
            # PBC on z
            if z == nb_vox[3] + 1
                z = 1
                y = y + 1
            end
            # PBC on y
            if y == nb_vox[2] + 1
                y = 1
                z = 1
                x = x + 1
            end
        end
    end

    # Updating volume
    volume = Volume( matrix, cell_matrix, nb_vox, center )

    # Scaling cell vectors by number of voxels
    # - Loop over dimensions
    for i=1:3
        cell_matrix[ :, i ] = cell_matrix[ :, i ]*nb_vox[i]
    end

    # Returns
    # - atoms_list: AtomList containing atomic positions, names
    # - cell_param: Cell_param contains cell information
    # - nb_vox: number of voxels in each directions
    # - matrix: data contained in the *.cube file
    # - vec_vox: voxel describing the voxel structure
    # - center: shift between position origin and grid origin
    return atom_list, matrix2Params( cell_matrix ), nb_vox, matrix, vec_vox, nb_vox, center
end
function readCubeVolume( file_path::T1 ) where { T1 <: AbstractString }
    # Argument
    # - file_path: path to the input file
    # Output
    # - atom_list: AtomList object with atom positions, names and index
    # - cell_matrix: cell matrix describing cell
    # - volume: volume with ELF, density, or whatever is contained in the cube file

    # Max number of column for data
    nb_col = 6

    # Open input file
    handle_in = open( file_path )

    # Reads entire files and put all lines in tuple
    lines = readlines( handle_in )

    # Closing input file
    close( handle_in )

    # Ignoring first 2 lines
    offset=3

    # Number of atoms
    nb_atoms = parse(Int, split( lines[ offset ] )[1] );

    # Origin position of the density
    center=zeros(3)
    for i=1:3
        center[i] = parse(Float64, split( lines[ offset ] )[ i + 1 ])*conversion.bohr2Ang
    end

    # Number of voxels in each direction
    nb_vox=zeros(Int,3)
    for i=1:3
        nb_vox[i] = parse(Float64, split( lines[ i + offset ] )[1] )
    end

    # Parse Cell information
    cell_matrix=zeros(Real, 3, 3 )
    for i=1:3
        for j=1:3
            cell_matrix[ i, j ] = parse(Float64, split( lines[ i + offset ] )[ j + 1 ] )*conversion.bohr2Ang
        end
    end

    # Increase offset
    offset = 6

    # Parse Atomic Information
    atom_list=atom_mod.AtomList(nb_atoms);

    # Loop over atoms
    for atom=1:nb_atoms
        # Get atom names
        atom_list.names[atom] = periodicTable.z2Names( parse(Int, split( lines[ atom + offset ] )[1] ) )

        # Loop over dimensions
        for j=1:3
            # Parsing
            atom_list.positions[ j, atom ] = parse(Float64, split( lines[ atom + offset ])[ j + 2 ] )*conversion.bohr2Ang
        end
    end

    # Parse Volume Data
    nb_tot = nb_vox[1]*nb_vox[2]*nb_vox[3]

    # Init matrix for data
    matrix = zeros(Real, nb_vox[1], nb_vox[2], nb_vox[3] )

    # Update offset
    offset = round(Int, offset + nb_atoms + 1 )

    # Init position of matrix
    x=1; y=1; z=1;

    # Loop over volume lines
    for i=0:nb_tot/nb_col-1
        # Loop over columns
        for j=1:nb_col
            # index of current line
            index_ = round(Int, offset + i  )
            # volume data
            keys = split( lines[ index_ ] )
            #
            if size( keys )[1] < j
                break
            end
            # Storing data in matrix
            matrix[ x, y, z ] = parse(Float64, keys[ j ] )
            # Increment z
            z = z + 1
            # PBC on z
            if z == nb_vox[3] + 1
                z = 1
                y = y + 1
            end
            # PBC on y
            if y == nb_vox[2] + 1
                y = 1
                z = 1
                x = x + 1
            end
        end
    end

    # Updating volume
    volume = Volume( matrix, cell_matrix, nb_vox, center )

    # Scaling cell vectors by number of voxels
    # - Loop over dimensions
    for i=1:3
        cell_matrix[ :, i ] = cell_matrix[ :, i ]*nb_vox[i]
    end

    # Returns atom, cell info and volume data
    return atom_list, matrix2Params( cell_matrix ), volume
end
#-----------------------------------------------------------------------------------

# Write information of a volume into a *.cube file
#-----------------------------------------------------------------------------------
function writeCube( file_path::T1, atoms::T2, cell::T3, volume::T4, title_line::T5="TITLE OF FRAME", comment_line::T6="COMMENT LINE OF FRAME" ) where { T1 <: AbstractString, T2 <: atom_mod.AtomList, T3 <: cell_mod.Cell_param, T4 <: cube_mod.Volume, T5 <: AbstractString, T6 <: AbstractString }
    # Arguments
    # - file_path: path to the output file
    # - atoms: AtomList contains atomic information
    # - cell: Cell_Param containing cell information
    # - volume: data from
    # - title_line (opt): Title of the frame
    # - comment_line (opt): Comment line describing cell
    # Output
    # - Bool: whether or not writting was successful

    # Maximum number of column for the volume data
    max_col = 6

    # Opening file
    handle_out = open( file_path, "w" )

    # Write title line
    Base.write( handle_out, string( title_line,   "\n" ) )

    # Write comment line
    Base.write( handle_out, string( comment_line, "\n" ) )

    # Write number of atoms
    Base.write( handle_out, string( size(atoms.positions)[2], " " ) )

    # Write origin vector
    # - Loop over dimension
    for i=1:3
        # Write vector component, rounding to 3 digits
        Base.write( handle_out, string( round( volume.origin[i]*conversion.ang2Bohr, digits=3 ), " " ) )
    end
    Base.write( handle_out, "\n" )

    # Writting Cell information
    # - Loop over dimensions
    for i=1:3
        # Write number of voxel in direction i
        Base.write( handle_out, string( volume.nb_vox[i], " " ) )
        # - Second loop over dimensions
        for j=1:3
            Base.write( handle_out, string( round( volume.vox_vec[i,j]*conversion.ang2Bohr/volume.nb_vox[i], digits=3 ), " " ) )
        end
        # Write end of line
        Base.write( handle_out, string("\n") )
    end

    # Writting atomic information
    # - Loop over atoms
    for atom=1:size( atoms.positions )[2]
        # Write atom name
        Base.write( handle_out, string( periodicTable.names2Z( atoms.names[atom] ), " " ) )
        # Write atom Z, again for some reason, with 5 digits if possible (to be fixed, but still works for VMD)
        Base.write( handle_out, string( round( periodicTable.names2Z( atoms.names[atom] ), digits=5 ), " " ) )
        # Loop over dimensions
        for i=1:3
            # Writting atomic positions (converting to atomic unit first)
            Base.write( handle_out, string( round( atoms.positions[i,atom]*conversion.ang2Bohr, digits=3 ), " " ) )
        end
        # Write end of line
        Base.write( handle_out, string("\n") )
    end

    # Writting volume info
    # - Loop over voxel in i dimension
    for i=1:volume.nb_vox[1]
        # - Loop over voxel in j dimension
        for j=1:volume.nb_vox[2]
            # - Loop over voxel in k dimension
            for k=1:volume.nb_vox[3]
                # Write volume data with rounding 5
                Base.write( handle_out, string( round( volume.matrix[i,j,k], digits=5 ), " " ) )
                # Make sure that we make maximum 6 columns
                if k % max_col == 0
                    Base.write( handle_out, string("\n") )
                end
            end
            # Write end of line
            write( handle_out, string("\n") )
        end
    end

    # Closing file
    close( handle_out )

    # Returns bool
    return true
end
#-----------------------------------------------------------------------------------

# Apply PBC on cube
#-----------------------------------------------------------------------------------
function applyPBCCube( index::T1, max_vox::T2 ) where { T1 <: Int, T2 <: Int }
    # Argument
    # - index: original index
    # - max_vox: maximum number of voxel in given direction
    # Output
    # - index: new index with PBC applied

    # Apply Lower PBC
    if index < 1
        index = max_vox + index
    end

    # Apply Upper PBC
    if index > max_vox
        index = index - max_vox
    end

    # Returns the modified index
    return index
end
#-----------------------------------------------------------------------------------

# Displace all values by a given vector
#-----------------------------------------------------------------------------------
function moveValues( volume::T1, vector_move::Vector{T2} ) where { T1 <: Volume, T2 <: Int }
    # Argument
    # - volume: Volume with all info
    # - vector_move: vector by which to move the data
    # Output
    # - volume: volume with shifted data

    # Initialize matrices
    matrix_in = copy( volume.matrix )
    matrix_out = copy( volume.matrix )

    # Loop over dimension 1
    for i=1:volume.nb_vox[1]
        # Handle very large move vectors
        if vector_move[1] > volume.nb_vox[1] || vector_move[1] < -volume.nb_vox[1]
            vector_move[1] = vector_move % volume.nb_vox[1]
        end

        # Computing move vector on direction i
        i_move = applyPBCCube( i + vector_move[1], volume.nb_vox[1] )

        # Loop over dimension 2
        for j=1:volume.nb_vox[2]
            # Handle very large move vectors
            if vector_move[2] > volume.nb_vox[2] || vector_move[2] < -volume.nb_vox[2]
                vector_move[2] = vector_move % volume.nb_vox[2]
            end

            # Computing move vector on direction i
            j_move = applyPBCCube( j + vector_move[2], volume.nb_vox[2] )

            # Loop over dimension 3
            for k=1:volume.nb_vox[3]
                # Handle very large move vectors
                if vector_move[3] > volume.nb_vox[3] || vector_move[3] < -volume.nb_vox[3]
                    vector_move[3] = vector_move % volume.nb_vox[3]
                end

                # Computing move vector on direction i
                k_move = applyPBCCube( k + vector_move[3], volume.nb_vox[3] )

                # Shifting matrix elements
                matrix_out[ i_move, j_move, k_move ] = matrix_in[ i, j, k ]
            end
        end
    end

    # Replacing matrix with new one
    volume.matrix = matrix_out

    # Return new volume with modified matrix
    return volume
end
#-----------------------------------------------------------------------------------

# Compute distance between two points in the grid, orthorombic cells only
#-----------------------------------------------------------------------------------
function distanceSimpleCube( position_1::Vector{T1}, position_2::Vector{T2}, volume::T3 ) where {T1 <: Int, T2 <: Int, T3 <: Volume }
    # Argument
    # - volume: information about the volume
    # - position_cube: position in the cube
    # - position_atom: position of the atom
    # Output
    # - dist: distance between atom and cube

    # Initialize distance
    dist = 0

    # Loop over dimensions
    for i=1:3
        # Compute local distance in each dimension
        dist_loc = ( position_1[i] - position_2[i] )*(norm( volume.vox_vec[i,:] )/volume.nb_vox[i] )
        # Square dimension distance
        dist = dist + dist_loc*dist_loc
    end

    # Return distance
    return sqrt(dist)
end
#-----------------------------------------------------------------------------------

# Keep only values within a given radius of positions
#-----------------------------------------------------------------------------
function carveCube( volume::T1, positions::Array{T2,2}, cut_off::T3, fac::T4=1, soft_fac::T5=0.5 ) where { T1 <: Volume, T2 <: Int, T3 <: Real, T4 <: Real, T5 <: Real }
    # Argument
    # - volume: volume information with data
    # - positions: positions of the atoms around which to carve the values
    # Output
    # - volume: modified volume with carved values

    # Check that the softening cut-off is >= 1 or set it to 1
    if fac < 1
        fac = 1.00005
    end

    # Check that the softening factor is between 0 and 1 or set it to 0.5 as default
    if soft_fac > 1 || soft_fac < 0
        soft_fac = 0.5
    end

    # Initialization of matrix
    matrix_in  = copy( volume.matrix )
    matrix_out = zeros( volume.nb_vox[1], volume.nb_vox[2], volume.nb_vox[3] )

    # Loop over dimension 1
    for i=1:volume.nb_vox[1]
        # Loop over dimension 2
        for j=1:volume.nb_vox[2]
            # Loop over dimension 3
            for k=1:volume.nb_vox[3]
                # Initialize bool to check whether to carve
                found  = false
                found2 = false

                # Loop over atoms
                for atom=1:size(positions)[2]
                    # Computes distance between positions on grid and atom position
                    dist = distanceSimpleCube( [ i, j, k ], positions[ :, atom ], volume )

                    # Check wether point is within a given radius
                    if  dist < cut_off
                        # If so keep matrix data
                        matrix_out[i,j,k] = matrix_in[i,j,k]
                        # Indicate that it's ok
                        found = true
                        break
                    # If the point is within a secondary cut-off we can damp the values
                    elseif dist < fac*cut_off
                        found2 = true
                    end
                end

                # Soften the carving so that the data looks less messy
                if ( ! found ) && found2
                    # Copy the data from original matrix, with a softening factor
                    # that is determined by user
                    matrix_out[i,j,k] = matrix_in[i,j,k]*soft_fac
                end
            end
        end
    end

    # Update volume with carved volume
    volume.matrix = matrix_out

    # Modified volume
    return volume
end
#-----------------------------------------------------------------------------

# Get index nearest to points
#-----------------------------------------------------------------------------
# Get the closest voxel node in the grid for a given atomic position
function getClosestIndex( position::Vector{T1}, volume::T2 ) where { T1 <: Real, T2 <: Volume }
    # Arguments
    # - position: position of the target atom
    # - volume: volume with all information
    # - cell: cell_param with all information
    # Output
    # - index: index closest of the target atom

    # Trying to guess the closest g
    index=zeros(Int,3)

    # Loop over the dimensions
    for i=1:3
        # Compute index in the ith dimension
        index[i] = round(Int, position[i]/(norm( volume.vox_vec[:,i] )/volume.nb_vox[i] ) ) + 1

        # Wrapping up index
        if index[i] > volume.nb_vox[i] || index[i] < - volume.nb_vox[i]
            index[i] = index[i] % volume.nb_vox[i]
        end

        # Applying PBC on index
        if index[i] < 1
            index[i] = volume.nb_vox[i] + index[i]
        end
    end

    # Returns the index
    return index
end
# Compute shift in voxels between origin of the voxel and position grids
function computeDisplacementOrigin( volume::T1 , cell::T2 ) where { T1 <: Volume, T2 <: cell_mod.Cell_param }
    # Argument
    # - volume: Volume containing all data
    # - cell: Cell_param with information about the cell
    # Output
    # - index: displacement in number of voxels of the origin

    # Initialize output vector
    index = zeros(Int, 3 )

    # Loop over dimension
    for i=1:3
        # Compute index in real
        index[i] = round(Int, volume.origin[i]/cell.length[i]*volume.nb_vox[i] )
    end

    # Returns the index of the vector
    return index
end
# Compute shift in voxels between origin of the voxel and position grids
function getClosestIndex( position::Vector{T1}, volume::T2 , cell::T3, origin_index::Vector{T4} ) where { T1 <: Real, T2 <: Volume, T3 <: cell_mod.Cell_param , T4 <: Int }
    # Argument
    # - position: position of atom to get into voxels
    # - volume: structure containing all data
    # - origin_index: shift of the voxel grid compared to position grid
    # Output
    # - index: position in voxel from the position in the cell

    # Trying to guess the closest grid point to the center
    index = zeros(Int, 3 )

    # Loop over dimensions
    for i=1:3
        # Compute original index
        index[i] = round(Int, position[i]/cell.length[i]*volume.nb_vox[i]) - origin_index[i]

        # Apply Upper PBC
        if index[i] > volume.nb_vox[i]
            index[i] = index[i] - volume.nb_vox[i]
        end

        # Apply lower PBC
        if index[i] < 1
            index[i] = volume.nb_vox[i]+index[i]
        end
    end

    # Return position in the voxel grid
    return index
end
# Get voxel point in the middle of two positions
function dataInTheMiddleWME( atoms::T1, cell::T2 , atom1::T3, atom2::T4, volume::T5 ) where { T1 <: atom_mod.AtomList, T2 <: cell_mod.Cell_param, T3 <: Int, T4 <: Int, T5 <: Volume }
    # Arguments
    # - atoms: AtomList containing all atomic information, including names, index and positions
    # - cell: Cell_param with all cell information
    # - atom1: index of first atom
    # - atom2: index of second atom
    # - volume: structure containing all informations
    # Output
    # - value of the volume at the center of atom 1 and atom 2

    # Copies positions of atom 1 and 2
    position1 = atoms.positions[ atom1,:]
    position2 = atoms.positions[ atom2,:]

    # Initialize data for positions between 1 and 2
    center = zeros(Real, 3 )

    # Moving 2 to closest image to 1 (can be out of the box)
    # - Loop over dimensions
    for i=1:3
        # - Compute differences in position between 1 and 2
        di = position1[i] - position2[i]

        # - Apply Upper PBC
        if di > cell.length[i]*0.5
            position2[i] = position2[i] + cell.length[i]
        end

        #  - Apply Lower PBC
        if di < -cell.length[i]*0.5
            position2[i] = position2[i] - cell.length[i]
        end

        # - Compute the position of the center (can be out of the box)
        center[i] = 0.5*(position1[i]+position2[i])

        # - Wrap center positions
        center[i] = cell_mod.wrap( center[i], cell.length[i] )
    end

    # Get indexes of the voxel of the center between 1 and 2
    index = getClosestIndex( center , data , cell )

    # Return value at the center
    return volume.matrix[ index[1], index[2], index[3] ]
end
# Get all values in the line between deux atomic positions
function traceLine( atom1::T1, atom2::T2, nb_points::T3, volume::T4, atoms::T5 , cell::T6 ) where { T1 <: Int, T2 <: Int, T3 <: Int, T4 <: Volume , T5 <: atom_mod.AtomList, T6 <: cell_mod.Cell_param }
    # Argument
    # - atom1, atom2: index of atoms 1 and 2
    # - nb_points: number of points to take in the line between 1 and 2
    # - volume: Volume structure containing volumic data
    # - atoms: AtomList containing atomic informations
    # - cell: Cell_Param containing all cell information
    # Output
    # - distances: distance from point 1
    # - data: data along the line from 1 to 2

    # Copy positions
    position1 = atoms.positions[atom1,:]
    position2 = atoms.positions[atom2,:]

    # Initialize positions from 1
    dp = zeros(Real, 3 )

    # Loop over dimensions
    for i=1:3
        # Wrapping positions
        position1[i] = cell_mod.wrap( position1[i], cell.length[i] )
        position2[i] = cell_mod.wrap( position2[i], cell.length[i] )

        # Moving 2 to closest image to 1 (can be out of the box)
        di = position1[i] - position2[i]

        # Apply Upper PBC
        if di > cell.length[i]*0.5
            position2[i] = position2[i] + cell.length[i]
        end

        # Apply lower PBC
        if di < -cell.length[i]*0.5
            position2[i] = position2[i] - cell.length[i]
        end

        # Get positions of point i
        dp[i] = ( position2[i] - position1[i] )/nb_points
    end

    dp_value = cell_mod.distance( atoms, cell, atom1, atom2 )/nb_points

    # Displacement due to origin of the volume
    origin=computeDisplacementOrigin( volume , cell )

    # Output Tables
    distances = zeros(Real, nb_points )
    data      = zeros(Real, nb_points )

    # Move cursor along the line and get data at each positions
    # - Start at position 1
    curseur=position1

    # - Loop over points
    for i=1:nb_points
        # Get the index of the current point
        indexs = getClosestIndex( curseur, volume, cell, origin )

        # Get data at current point
        data[i] = volume.matrix[ indexs[1], indexs[2], indexs[3] ]

        # Move cursor
        for j=1:3
            curseur[j] += dp[j]
        end

        # Compute the distance from 1
        if i > 1
            # Loop over dimensions
            for j=1:3
                # Add distance along all point
                distances[i] = distances[i-1] + dp_value
            end
        end
    end

    # Return distances from 1 to 2 and associated data along the line
    return distances, data
end
#-----------------------------------------------------------------------------

end
