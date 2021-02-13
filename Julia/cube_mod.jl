module cube_mod

using atom_mod
using cell_mod
using geom

export readCube, dataInTheMiddleWME, getClosestIndex, traceLine

#-----------------------------------------------------------------------
# Volume 4, nb_vox*nb_vox*nb_vox
#-----------------------
# 1-3 for positions
# 4 for the values
#-----------------------
mutable struct Volume
    matrix::Array{Real}
    vox_vec::Array{Real}
    nb_vox::Vector{Int}
    origin::Vector{Real}
    function Volume()
        new( Array{ Real, 3 }(undef,1,1,1), Array{Real,2}(undef,3,3), Array{Int,1}(undef,3),Array{Real,1}(undef,3) )
    end
    function Volume( nb_vox_iso::T1 ) where {T1 <: Real}
        if nb_vox_iso > 0
            new( Array{Real}( nb_vox_iso, nb_vox_iso, nb_vox_iso), Array{Real}(3,3), [nb_vox_iso, nb_vox_iso, nb_vox_iso] )
        else
            print("/!\\ Error: Wrong size for the number of voxels.");
        end
    end
    function Volume( nb_vox::Vector{T1} ) where {T1 <: Real}
        if size(nb_vox)[1] == 3
            new( Array{Real}(nb_vox[1],nb_vox[2],nb_vox[3]),Array{Real}(3,3),nb_vox);
        else
            print("/!\\ Error: Wrong size for the number of voxels.");
        end
    end
end
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------------------
# Reads a cube file and returns all or parts of its informations
function readCube( file_name::T1 ) where { T1 <: AbstractString }
    #--------------
    # Reading file
    #----------------------
    file=open(file_name);
    lines=readlines(file);
    close(file);
    #-----------------------

    #------------------
    # Number of atoms
    #-----------------------------------------
    nb_atoms = parse(Int,split(lines[3])[1]);
    #-----------------------------------------

    #--------------------------------
    # Origin position of the density
    #-----------------------------------------------------
    center=zeros(3)
    for i=1:3
        center[i]=parse(Float64,split(lines[3])[1+i])*0.529177
    end
    #-----------------------------------------------------

    #-------------------------------------
    # Number of voxels in each direction
    #-----------------------------------------------------
    nb_vox=zeros(Int,3)
    for i=1:3
        nb_vox[i] = parse(Float64, split( lines[3+i] )[1] )
    end
    #-----------------------------------------------------

    #--------------------
    # Reads Cell Matrix
    #-----------------------------------------------------
    cell_matrix=zeros(Real,3,3)
    for i=1:3
        for j=1:3
            cell_matrix[i,j] = parse(Float64, split( lines[3+i] )[1+j] )*0.529177
        end
    end
    #-----------------------------------------------------

    #-------------
    # Reads atoms
    #----------------------------------------------------
    atom_list=atom_mod.AtomList(nb_atoms);
    for i=1:nb_atoms
        atom_list.names[i] = split( lines[6+i] )[1]
        for j=1:3
            atom_list.positions[i,j] = parse(Float64, split( lines[6+i])[2+j] )*0.529177
        end
    end
    #----------------------------------------------------

    #-----------------------
    # Reads Volumetric Data
    #----------------------------------------------------
    nb_tot=nb_vox[1]*nb_vox[2]*nb_vox[3]
    nb_col=6
    matrix = Array{Real,3}( undef,nb_vox[1], nb_vox[2], nb_vox[3] )
    offset=(Int)(6+nb_atoms+1)
    x=1; y=1; z=1;
    for i=0:nb_tot/nb_col-1
        for j=1:nb_col
            matrix[x,y,z] = parse(Float64, split( lines[(Int)(offset+i)])[j] )
            z=z+1;
            if z == nb_vox[3]+1
                z=1;
                y=y+1;
            end
            if y == nb_vox[2]+1
                y=1;
                z=1;
                x=x+1;
            end
        end
    end
    #----------------------------------------------------

    #-----------------
    # Updating volume
    #--------------------------------
    volume=Volume()
    volume.nb_vox=nb_vox
    volume.vox_vec=cell_matrix
    volume.matrix=matrix
    volume.origin=center
    #---------------------------------

    #-------------------------------------------
    # Scaling cell vectors by number of voxels
    #--------------------------------------------
    cell_vecs=cell_mod.Cell_matrix(cell_matrix)
    for i=1:3
        cell_vecs.matrix[i,i]=cell_vecs.matrix[i,i]*nb_vox[i]
    end
    #---------------------------------------------

    return atom_list, cell_matrix, volume
end
#-----------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------
function computeDisplacementOrigin( data::T1 , cell::T2 ) where { T1 <: Volume, T2 <: cell_mod.Cell_param }
    index=zeros(Int,3)
    for i=1:3
        guess=data.origin[i]/cell.length[i]*data.nb_vox[i]
        index[i]=trunc(guess)
        if guess-index[i] > 0.5
            index[i] += 1
        end
    end
    return index
end

function getClosestIndex( position::Vector{T1}, volume::T2 , cell::T3 ) where { T1 <: Real, T2 <: Volume, T3 <: cell_mod.Cell_param }
    # Trying to guess the closest grid point to the center
    index=zeros(Int,3)
    # Rounding up the raw guess
    for i=1:3
        index[i] = trunc(position[i]/cell.length[i]*volume.nb_vox[i])  + 1
    end
    # Displacement due to origin
    mod=computeDisplacementOrigin( volume , cell )
    for i=1:3
        index[i]-=mod[i]
    end
    for i=1:3
        if index[i] > volume.nb_vox[i]
            index[i] = index[i] - volume.nb_vox[i]
        end
        if index[i] < 1
            index[i] = volume.nb_vox[i]+index[i]
        end
    end
    return index
end

function getClosestIndex( position::Vector{T1}, volume::T2 , cell::T3, origin_index::Vector{T4} ) where { T1 <: Real, T2 <: Volume, T3 <: cell_mod.Cell_param , T4 <: Int }
    # Trying to guess the closest grid point to the center
    index=zeros(Int,3)
    # Rounding up the raw guess
    for i=1:3
        index[i] = trunc(position[i]/cell.length[i]*volume.nb_vox[i])  + 1
    end
    for i=1:3
        index[i]-= origin_index[i]
    end
    for i=1:3
        if index[i] > volume.nb_vox[i]
            index[i] = index[i] - volume.nb_vox[i]
        end
        if index[i] < 1
            index[i] = volume.nb_vox[i]+index[i]
        end
    end
    return index
end

function dataInTheMiddleWME( atoms::T1, cell::T2 , atom1::T3, atom2::T4, data::T5 ) where { T1 <: atom_mod.AtomList, T2 <: cell_mod.Cell_param, T3 <: Int, T4 <: Int, T5 <: Volume }
    # Wrapped  Edition
    # Copies of 1 and 2
    position1 = atoms.positions[atom1,:]
    position2 = atoms.positions[atom2,:]
    # Moving 2 to closest image to 1 (can be out of the box)
    for i=1:3
        di = position1[i] - position2[i]
        if di > cell.length[i]*0.5
            position2[i] = position2[i] + cell.length[i]
        end
        if di < -cell.length[i]*0.5
            position2[i] = position2[i] - cell.length[i]
        end
    end
    # compute the position of the center (can be out of the box)
    center=zeros(Real,3)
    for i=1:3
        center[i] = 0.5*(position1[i]+position2[i])
    end
    # wrap the center
    for i=1:3
        center[i] = cell_mod.wrap( center[i], cell.length[i] )
    end
    index=getClosestIndex( center , data , cell )
    return data.matrix[index[1],index[2],index[3]]
end

# Trace the volume between two points.
function traceLine( atom1::T1, atom2::T2, nb_points::T3, volume::T4, atoms::T5 , cell::T6 ) where { T1 <: Int, T2 <: Int, T3 <: Int, T4 <: Volume , T5 <: atom_mod.AtomList, T6 <: cell_mod.Cell_param }

    # Extracting positions
    position1 = atoms.positions[atom1,:]
    position2 = atoms.positions[atom2,:]

    # Wrapping
    for i=1:3
        position1[i]=cell_mod.wrap(position1[i],cell.length[i])
        position2[i]=cell_mod.wrap(position2[i],cell.length[i])
    end

    # Moving 2 to closest image to 1 (can be out of the box)
    for i=1:3
        di = position1[i] - position2[i]
        if di > cell.length[i]*0.5
            position2[i] = position2[i] + cell.length[i]
        end
        if di < -cell.length[i]*0.5
            position2[i] = position2[i] - cell.length[i]
        end
    end

    # Move vector and distance
    dp=zeros(Real,3)
    for i=1:3
        dp[i]=(position2[i]-position1[i])/nb_points
    end
    dp_value=cell_mod.distance(atoms,cell,atom1,atom2)/nb_points

    # Displacement due to origin of the volume
    origin=computeDisplacementOrigin( volume , cell )

    # Output Tables
    distances=zeros(Real,nb_points)
    data=zeros(Real,nb_points)

    # Moving along the lines
    curseur=position1
    for i=1:nb_points
        indexs=getClosestIndex(curseur,volume,cell,origin)
        # Data
        data[i]=volume.matrix[indexs[1],indexs[2],indexs[3]]
        # Movement
        for j=1:3
            curseur[j] += dp[j]
        end
        if i > 1
            for j=1:3
                distances[i]=distances[i-1]+dp_value
            end
        end
    end

    return distances, data
end

end
