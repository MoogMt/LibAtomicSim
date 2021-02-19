module exp_data

# Load necessary modules from LibAtomicSim
using atom_mod
using cell_mod
using filexyz
using pdb
using fftw
using correlation
using conversion

# Description
# - Functions used to compute experimental data from simulations

# Exporting functions
export vdosFromPosition
export computeGr
export computeFQ
export computeBarycenter
export computeMSD
export computeRMSD

# TODO:
# Recheck everything for consistency

# Vibration Density of States
#==============================================================================#
# file_traj: path to the TRAJEC.xyz file
# max_lag_frac: max fraction of total time to be used for correlation function (0<x<0.5)
# dt: timestep in ps (important)
function vdosFromPosition( file_traj::T1 , max_lag_frac::T2 , dt::T3 ) where { T1 <: AbstractString, T2 <: Real, T3 <: Real }

        # Reading Trajectory
        traj = filexyz.readFileAtomList(file_traj)
        if traj == false
            return false, false
        end

        # Computing velocities
        dx=1
        velocity=cell_mod.velocityFromPosition(traj,dt,dx)

        nb_atoms=size(traj[1].names)[1]
        nb_step=size(velocity)[1]

        # Compute scalar product
        velo_scal=zeros(nb_step,nb_atoms)
        for atom=1:nb_atoms
            for step=1:nb_step
                for i=1:3
                    velo_scal[step,atom] += velocity[step,atom,i]*velocity[1,atom,i]
                end
            end
        end

        max_lag=Int(trunc(nb_step*max_lag_frac))

        # Average correlation
        freq=zeros(trunc(Int,max_lag/2))
        vdos=zeros(trunc(Int,max_lag/2))
        for atom=1:nb_atoms
            freq, vdos_loc = fftw.doFourierTransformShift( correlation.autocorrNorm( velo_scal[:,atom] , max_lag ), dt )
            for i=1:size(vdos_loc)[1]
                vdos[i] += vdos_loc[i]
            end
        end
        for i=1:size(vdos)[1]
            vdos[i] = vdos[i]/nb_atoms
        end
        # Conversion to cm-1
        for i=1:size(freq)[1]
            freq[i]=freq[i]*conversion.tHz2cm
        end

    return freq, vdos
end
# Same thing but write the results in file_out
function vdosFromPosition( file_traj::T1 , file_out::T2 , max_lag_frac::T3 , dt::T4 ) where { T1 <: AbstractString, T2 <: AbstractString, T3 <: Real, T4 <: Real }

    freq, vdos = vdosFromPosition( file_traj , max_lag_frac , dt )
    if freq == false
        return false, false
    end

    # Writting data to file
    file_o=open(string(file_out),"w")
    for i=1:size(vdos)[1]
        Base.write(file_o,string(freq[i]," ",vdos[i],"\n"))
    end
    close(file_o)

    return freq, vdos
end
# Same thing but divides trajectory in nb_windows time windows
# and averages the resulting VDOS for smoother result
function vdosFromPosition( file_traj::T1 , max_lag_frac::T2 , dt::T3, nb_windows::T4 ) where { T1 <: AbstractString, T2 <: Real, T3 <: Real, T4 <: Int }

        # Reading Trajectory
        traj = filexyz.readFileAtomList(file_traj)
        if traj == false
            return false, false
        end

        # Computing velocities
        dx=1
        velocity=cell_mod.velocityFromPosition(traj,dt,dx)

        nb_atoms=size(velocity)[2]
        nb_step=Int(trunc(size(velocity)[1]/nb_windows))

        # Compute scalar product
        velo_scal=zeros(nb_step,nb_atoms,nb_windows)
        for atom=1:nb_atoms
            for window=1:nb_windows
                start_win=(window-1)*nb_step+1
                end_win=window*nb_step-1
                count_=1
                for step=start_win:end_win
                    for i=1:3
                        velo_scal[count_,atom,window] += velocity[step,atom,i]*velocity[1,atom,i]
                    end
                    count_+=1
                end
            end
        end
        max_lag=Int(trunc(nb_step*max_lag_frac))

        # Average correlation
        freq=zeros(Int(trunc(max_lag*max_lag_frac))-1)
        vdos=zeros(Int(trunc(max_lag*max_lag_frac))-1)
        for window=1:nb_windows
            for atom=1:nb_atoms
                freq, vdos_loc = fftw.doFourierTransformShift( correlation.autocorrNorm( velo_scal[:,atom,window] , max_lag ), dt )
                for i=1:size(vdos)[1]
                    vdos[i] += vdos_loc[i]
                end
            end
            vdos /= nb_atoms
        end

        # Conversion to cm-1
        freq=freq.*conversion.tHz2cm

    return freq, vdos
end
# Same thing but write the results in file_out
function vdosFromPosition( file_traj::T1 , file_out::T2 , max_lag_frac::T3 , dt::T4, nb_windows::T5 ) where { T1 <: AbstractString, T2 <: AbstractString, T3 <: Real, T4 <: Real, T5 <: Int }

    freq, vdos = vdosFromPosition( file_traj , max_lag_frac , dt, nb_windows )
    if freq == false
        return false, false
    end

    # Writting data to file
    file_o=open(string(file_out),"w")
    for i=1:size(vdos)[1]
        Base.write(file_o,string(freq[i]," ",vdos[i],"\n"))
    end
    close(file_o)

    return freq, vdos
end
#==============================================================================#

# Radial Distribution Function
#==============================================================================#
function computeGr( file_in::T1, a::T2, rmin::T3, rmax::T4, dr::T5 ) where { T1 <: AbstractString, T2 <: Real, T3 <: Real, T4 <: Real, T5 <: Real }

    # Reading trajectory
    traj = filexyz.readFileAtomList(file_in)
    if traj == false
        return false
    end

    #Building cell
    cell=cell_mod.Cell_param(a,a,a)

    # Volume
    V=a^3

    # Getting nb atoms and steps
    nb_atoms=size(traj[1].names)[1]
    nb_step=size(traj)[1]
    nb_box=round(Int,(rmax-rmin)/dr)+1
    gr=zeros(nb_box)

    # Computing g_r
    for step=1:nb_step
        print("G(r) Progress: ",step/nb_step*100,"\n")
        for atom1=1:nb_atoms
            for atom2=atom1+1:nb_atoms
                dist=cell_mod.distance(traj[step],cell,atom1,atom2)
                if dist < rmax
                    gr[ Int( trunc( (dist-rmin)/dr ) )+1  ] += 1
                end
            end
        end
    end
    # Normalization
    for i=1:nb_box
        gr[i] = gr[i]*V/(nb_step*nb_atoms*(nb_atoms-1)/2)/(4*pi*dr*(i*dr+rmin)^2)
    end

    return gr
end
function computeGr( file_in::T1, file_out::T2, V::T3, rmin::T4, rmax::T5, dr::T6 ) where { T1 <: AbstractString, T2 <: AbstractString, T3 <: Real, T4 <: Real, T5 <: Real, T6 <: Real }

    gr = computeGr(file_in,V,rmin,rmax,dr)

    if gr == false
        return false
    end

    file_o=open(file_out,"w")
    for i=1:size(gr)[1]
        Base.write(file_o,string(i*dr+rmin," ",gr[i],"\n"))
    end
    close(file_o)

    return gr
end
#==============================================================================#

# Structure Factor
#==============================================================================#
# Do not trust this ( MM 03/2020 )
function computeFQ( gr::Vector{T1}, rmin::T2, rmax::T3, dr::T4, rho::T5 ) where { T1 <: Real, T2 <: Real, T3 <: Real, T4 <: Real, T5 <: Real }
    nb_box=1000
    nb_box_int=Int(trunc((rmax-rmin)/dr))
    fq=zeros(nb_box)
    q=zeros(nb_box)
    dq=0.01
    for i=1:nb_box
        loc=0
        q[i] = i*dq
        for j=1:nb_box_int-1
            r=j*dr+rmin
            r2=(j+1)*dr+rmin
            loc += (sin(q[i]*r)*r*(gr[j]-1)+sin(q[i]*r2)*r2*(gr[j+1]-1))*0.5*dr
        end
        loc *= 4*pi/q[i]*rho
        fq[i] = 1 + loc
    end
    return q, fq
end
function computeFQ( file_out::T1, gr::Vector{T2}, rmin::T3, rmax::T4, dr::T5, rho::T6 ) where { T1 <: AbstractString, T2 <: Real, T3 <: Real, T4 <: Real, T5 <: Real, T6 <: Real }
    q,fq=computeFQ(gr,rmin,rmax,dr,rho)

    file_o=open(file_out,"w")
    for i=1:size(fq)[1]
        if q[i] > 0
            Base.write(file_o,string(q[i]," ",fq[i],"\n"))
        end
    end
    close(file_o)

    return q,fq
end
#==============================================================================#

# Barycenter
#==============================================================================#
function computeBarycenter( positions::Array{T1,2} ) where { T1 <: Real }
    barycenter=zeros(3)
    nb_atoms=size(positions)[1]
    for i=1:3
        for atom=1:nb_atoms
            barycenter[i] += positions[atom,i]
        end
    end
    barycenter /= nb_atoms
    return barycenter
end
function computeBarycenter( positions::Array{T1,2}, index_types::Vector{T2} ) where { T1 <: Real, T2 <: Int }
    barycenter=zeros(3)
    nb_atoms=size(positions)[1]
    for i=1:3
        for atom in index_types
            barycenter[i] += positions[atom,i]
        end
    end
    barycenter /= nb_atoms
    return barycenter
end
function computeBarycenter( positions::Array{T1,3} ) where { T1 <: Real }
    nb_step=size(positions)[1]
    nb_atoms=size(positions)[2]
    barycenter=zeros(nb_step,3)
    for step=1:nb_step
        barycenter[step,:] = computeBarycenter( positions[step,:,:] )
    end
    barycenter /= nb_atoms
    return barycenter
end
function computeBarycenter( positions::Array{T1,3}, index_types::Vector{T2} ) where { T1 <: Real, T2 <: Int }
    nb_step=size(positions)[1]
    nb_atoms=size(positions)[2]
    barycenter=zeros(nb_step,3)
    for step=1:nb_step
        barycenter[step,:]=computeBarycenter( positions[step,:,:], index_types )
    end
    barycenter /= nb_atoms
    return barycenter
end
function computeBarycenter( positions::Array{T1,3}, types::Vector{T2}, types_names::Vector{T3}, type_masses::Vector{T4} ) where { T1 <: Real, T2 <: AbstractString, T3 <: AbstractString, T4 <: Real }
    nb_step=size(positions)[1]
    nb_atoms=size(positions)[2]
    nb_types=size(types_names)[1]
    barycenter_ = zeros(nb_step,3)
    total_mass=0
    for type_=1:nb_types
        index_types=atom_mod.getTypeIndex(types,types_names[type_])
        nb_atoms_type=size(index_types)[1]
        barycenter_ += type_masses[type_]*nb_atoms_type*computeBarycenter(positions,index_types)
        total_mass += type_masses[type_]*nb_atoms_type
    end
    barycenter_ /= (total_mass)
    return barycenter_
end
#==============================================================================#

# MSD
#==============================================================================#
function computeMSD( positions::Array{T1,2}, barycenter_global::Array{T2,2} ) where { T1 <: Real, T2 <: Real }
    nb_step=size(positions)[1]
    msd=zeros(nb_step)
    for step=1:nb_step
        dist=0
        for i=1:3
            dist_loc=( positions[step,i] - barycenter_global[ step, i ] ) - ( positions[ 1, i ] - barycenter_global[ 1, i ] )
            dist += dist_loc*dist_loc
        end
        msd[step] = dist
    end
    return msd
end
function computeMSD( positions::Array{T1,3}, barycenter_global::Array{T2,2} ) where { T1 <: Real, T2 <: Real }
    nb_step  = size( positions )[1]
    nb_atoms = size( positions )[2]
    msd = zeros( nb_step )
    for step=1:nb_step
        for atom=1:nb_atoms
            dist=0
            for i=1:3
                dist_loc = ( positions[ step, atom, i] - barycenter_global[ step, i ] ) - ( positions[ 1, atom, i ] - barycenter_global[ 1, i ] )
                dist += dist_loc*dist_loc
            end
            msd[step] += dist
        end
        msd[step] /= nb_atoms
    end
    return msd
end
function computeMSD( traj::Vector{T1}, names::Vector{T2}, masses::Vector{T3} ) where { T1 <: AtomList, T2 <: AbstractString, T3 <: Real }
    atom_name_list=traj[1].names
    nb_step=size(traj)[1]
    # Computing barycenter movement
    barycenter_global=computeBarycenter( atom_mod.getPositionsAsArray( traj ), atom_name_list, names, masses )
    for step=1:nb_step
        barycenter_global[step,:] = barycenter_global[step,:] - barycenter_global[1,:]
    end

    positions = atom_mod.getPositionsAsArray( traj )

    # Computing MSD
    msd=zeros(nb_step)
    nb_atoms=0
    nb_types=size(names)[1]
    for name in names
        index_types=atom_mod.getTypeIndex( atom_name_list, name )
        nb_atoms += size(index_types)[1]
        for atom in index_types
            msd += computeMSD( positions[:,atom,:], barycenter_global )
        end
    end

    return msd/nb_atoms
end
function computeMSD( traj::Vector{T1}, names::Vector{T2}, masses::Vector{T3},  nb_window::T4 ) where { T1 <: AtomList, T2 <: AbstractString, T3 <: Real, T4 <: Int }
    atom_name_list=traj[1].names
    nb_step=size(traj)[1]
    # Computing barycenter movement
    barycenter_global=computeBarycenter(atom_mod.getPositions(traj),atom_name_list,names,masses)
    for step=1:nb_step
        barycenter_global[step,:] = barycenter_global[step,:] - barycenter_global[1,:]
    end
    # Extract positions
    positions=atom_mod.getPositions(traj)
    # Computing MSD
    msd=zeros(nb_step)
    nb_atoms=0
    nb_types=size(names)[1]
    for name in names
        index_types=atom_mod.getTypeIndex( atom_name_list, name )
        nb_atoms += size(index_types)[1]
        for atom in index_types
            msd += computeMSD( positions[:,atom,:], barycenter_global )
        end
    end
    return msd/nb_atoms
end
#==============================================================================#

# RMSD
#==============================================================================#
function computeRMSD( structure1::T1, structure2::T2 ) where { T1 <: atom_mod.AtomList, T2 <: AtomList }
    nb_atoms1=size(structure1.names)[1]
    nb_atoms2=size(structure2.names)[1]
    if nb_atoms1 != nb_atoms2
        print("Comparing two structures with different number of atoms. Stopping now!")
        return -1
    end
    rmsd=0
    for atom=1:nb_atoms1
        for i=1:3
            dist=structure1.positions[atom,i]-structure2.positions[atom,i]
            rmsd += dist*dist
        end
    end
    return rmsd/nb_atoms1
end
#==============================================================================#

end
