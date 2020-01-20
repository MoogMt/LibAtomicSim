module markov

using atom_mod
using cell_mod
using cube_mod
using clustering
using contact_matrix

export buildCoordinationMatrix, assignDataToStates, isolateSignificantStates
export transitionMatrix, chappmanKormologov
export writeStates, writeStateMatrix, writeTransitionsMatrix
export readStates, readTransitionMatrix

function buildCoordinationMatrix( traj::Vector{T1}, cell::T2, cut_off_bond::T3, max_neighbour::T4 ) where { T1 <: atom_mod.AtomList, T2 <: cell_mod.Cell_param , T3 <: Real, T4 <: Int }
    # general information about the simulation
    nb_atoms=size(traj[1].names)[1]
    nb_steps=size(traj)[1]
    # Type stuff
    types=[traj[1].names[1]]
    types_number=ones(1)
    count_types=zeros(Int,1)
    for i=1:nb_atoms
        found=false
        for j=1:size(types)[1]
            if types[j] == traj[1].names[i]
                count_types[j] += 1
                found=true
            end
        end
        if ! found
            push!(types,traj[1].names[i])
            push!(types_number,size(types)[1])
            push!(count_types,1)
        end
    end
    nb_type=size(types)[1]
    type_list=zeros(Int,nb_atoms)
    for i=1:nb_atoms
        for j=1:nb_type
            if traj[1].names[i] == types[j]
                type_list[i]=j
            end
        end
    end
    # Actual computation of the coordination matrix
    coord_matrix=zeros(nb_atoms,nb_steps,max_neighbour*nb_type)
    for step_sim=1:nb_steps
        print("Building Coordination Signal - Progress: ",step_sim/nb_steps*100,"%\n")
        # Bond Matrix
        bond_matrix=zeros(nb_atoms,nb_atoms)
        for atom1=1:nb_atoms
            for atom2=atom1+1:nb_atoms
                if cell_mod.distance( traj[step_sim], cell, atom1, atom2) < cut_off_bond
                    bond_matrix[atom1,atom2]=1
                    bond_matrix[atom2,atom1]=1
                end
            end
        end
        # Cleaning weird 3-member rings
        for atom1=1:nb_atoms-1
            for atom2=atom1+1:nb_atoms
                if bond_matrix[atom1,atom2] == 1
                    # For each bonded pair of atoms, check whether they are
                    # both bonded to another atom, forming unatural 3-member
                    # ring
                    for atom3=1:nb_atoms
                        if atom3 == atom1 || atom3 == atom1
                            continue
                        end
                        if bond_matrix[atom3,atom1] == 1 && bond_matrix[atom3,atom2] == 1
                            if cell_mod.distance(traj[step_sim],cell,atom1,atom3) > cell_mod.distance(traj[step_sim],cell,atom1,atom2) && cell_mod.distance(traj[step_sim],cell,atom1,atom3) > cell_mod.distance(traj[step_sim],cell,atom2,atom3)
                                bond_matrix[atom1,atom3]=0
                                bond_matrix[atom3,atom1]=0
                            elseif cell_mod.distance(traj[step_sim],cell,atom1,atom2) > cell_mod.distance(traj[step_sim],cell,atom2,atom3)
                                bond_matrix[atom1,atom2]=0
                                bond_matrix[atom2,atom1]=0
                            else
                                bond_matrix[atom2,atom3]=0
                                bond_matrix[atom3,atom2]=0
                            end
                        end
                    end
                end
            end
        end
        # Compute coord matrix
        for atom1=1:nb_atoms
            for type=1:nb_type
                count_coord=(type-1)*max_neighbour+1
                for atom2=1:nb_atoms
                    if (atom1 == atom2) || (types[type] != traj[step_sim].names[atom2])
                        continue
                    end
                    if bond_matrix[atom1,atom2] > 0
                        coord_matrix[atom1,step_sim,count_coord]=sum(bond_matrix[atom2,:])
                        count_coord += 1
                    end
                    if count_coord == type*max_neighbour
                        break
                    end
                end
                #sorting by type
                for i=(type-1)*max_neighbour+1:type*max_neighbour-1
                    #print("type = ",type,"; i = ",i)
                    for j=i+1:type*max_neighbour
                        #print(" j = ",j," ")
                        if coord_matrix[atom1,step_sim,i] < coord_matrix[atom1,step_sim,j]
                            stock = coord_matrix[atom1,step_sim,j]
                            coord_matrix[atom1,step_sim,j]=coord_matrix[atom1,step_sim,i]
                            coord_matrix[atom1,step_sim,i]=stock
                        end
                    end
                    #print("\n")
                end
            end
        end
    end
    return coord_matrix, types, type_list
end
function assignDataToStates( data::Array{T1,3}, nb_types::T2, types_number::Vector{T3} ) where { T1 <: Real, T2 <: Int , T3 <: Int }

    nb_series = size(data)[1]
    nb_steps  = size(data)[2]
    dim_data  = size(data)[3]

    number_per_types=zeros(Int,nb_types)
    for i=1:nb_types
        for j=1:nb_series
            if types_number[j] == i
                number_per_types[i] += 1
            end
        end
    end

    states=[zeros(Int,0,dim_data)]
    counts=[zeros(Real,0)]
    states_matrices=[zeros(Int,number_per_types[1],nb_steps)]
    for i=2:nb_types
        push!(states_matrices,zeros(number_per_types[i],nb_steps))
        push!(counts,zeros(Real,0))
        push!(states,zeros(Int,0,dim_data))
    end

    for type=1:nb_types
        for serie=1:number_per_types[type]
            for step=1:nb_steps
                print("State assignement - Type: ",type," - Series progress : ",serie/number_per_types[type]*100,"% - Step progress: ",step/nb_steps*100,"%\n")
                # Initiatilization
                nb_states=size(states[type])[1]
                if nb_states == 0
                    states_matrices[type][serie,step] = size(states[type])[1] + 1
                    states[type]=vcat(states[type],transpose(data[serie+sum(number_per_types[1:type-1]),step,:]))
                    push!(counts[type],1)
                else
                    found=false
                    # Loop over states
                    for state=1:nb_states
                        dist=0
                        for i=1:dim_data
                            dist += (states[type][state,i]-data[serie+sum(number_per_types[1:type-1]),step,i])*(states[type][state,i]-data[serie+sum(number_per_types[1:type-1]),step,i])
                        end
                        if dist == 0
                            counts[type][state] += 1
                            states_matrices[type][serie,step] = state
                            found=true
                            break
                        end
                    end
                    if ! found
                        states_matrices[type][serie,step] = size(states[type])[1] + 1
                        states[type]=vcat(states[type],transpose(data[serie+sum(number_per_types[1:type-1]),step,:]))
                        push!(counts[type],1)
                    end
                end
            end
        end
    end

    return states, states_matrices, counts
end
function assignDataToStates( data::Array{T1,3}, states::Array{T2,2} , nb_types::T3 , type_states::Vector{T4}, type_atoms::Vector{T5}, Err::T6 ) where { T1 <: Real, T2 <: Real , T3 <: Int, T4 <: Int, T5 <: Int, T6 <: Bool }
    nb_series = size(data)[1]
    nb_steps  = size(data)[2]
    dim_data  = size(data)[3]
    nb_states = size(states)[1]
    state_matrix=zeros(Int, nb_series, nb_steps )
    percent_states=zeros(Real, nb_states )
    unused=0

    for j=1:nb_steps
        print("Assigning data to states - Progress: ",j/nb_steps*100,"%\n")
        for i=1:nb_series
            for l=1:nb_states
                if type_atoms[i] != type_states[l]
                    continue
                end
                d=0
                for k=1:dim_data
                    d+= ( data[i,j,k] - states[l,k])*(data[i,j,k] - states[l,k])
                end
                if d == 0
                    percent_states[l] = percent_states[l] + 1
                    state_matrix[i,j] = l
                    break
                end
            end
            if state_matrix[i,j] == -1
                if Err
                    print("Out ",i," ",j," ",dim_data," ")
                    for k=1:dim_data
                        print(data[i,j,k]," ")
                    end
                    print("\n")
                end
                unused += 1
            end
        end
    end

    # Counting the number of participant states for each type
    type_count=zeros(nb_types)
    for type=1:nb_types
        for state=1:nb_states
            if type_states[state] == type
                type_count[type] += percent_states[state]
            end
        end
    end

    # Normalizing
    for state=1:nb_states
         percent_states[state] = percent_states[state]/type_count[type_states[state]]*100
    end

    return state_matrix, percent_states
end
function isolateSignificantStates( old_states::Array{T1,2}, percent_states::Vector{T2}, cut_off_states::T3, type_states::Vector{T4} ) where { T1 <: Real, T2 <: Real, T3 <: Real,  T4 <: Int }
    old_nb_states=size(old_states)[1]
    dim=size(old_states)[2]
    new_nb_states=0
    states_kept=zeros(0,dim)
    types_states_kept=zeros(Int,0)
    for i=1:old_nb_states
        if percent_states[i] > cut_off_states
            states_kept=[ states_kept ; transpose( old_states[i,:]) ]
            push!(types_states_kept,type_states[i])
        end
    end
    return states_kept, types_states_kept
end
function transitionMatrix( states::Array{T1,2}, state_matrix::Vector{T2}, nb_type::T4, type_series::Vector{T5}, min_lag::T6, max_lag::T7, d_lag::T8) where { T1 <: Int, T2 <: Real, T4 <: Real, T5 <: Real, T6 <: Real, T7 <: Int, T8 <: Int }

    nb_states=size(states)[1]
    nb_steps=size(state_matrix)[1]
    nb_lag_points=Int(trunc((max_lag-min_lag)/d_lag))

    state_transition_matrix=zeros( nb_states, nb_states, nb_lag_points )
    count_lag=1
    for lag=min_lag:d_lag:max_lag-1
        print("Computing Transition Matrix - Progress: ",lag/max_lag*100,"%\n")
        for j=lag+1:nb_steps
            if state_matrix[j-lag] == 0 || state_matrix[j] == 0
                continue
            end
            state_transition_matrix[ state_matrix[j-lag], state_matrix[j], count_lag ] += 1
        end
        count_lag += 1
    end

    # Normalization
    for lag=1:nb_lag_points
        print("Normalizing Transition Matrix: ",lag/nb_lag_points*100,"%\n")
        for i=1:nb_states
            sum_transition=sum( state_transition_matrix[i,:,lag] )
            if sum_transition != 0
                state_transition_matrix[i,:,lag] /= sum_transition
            end
        end
    end

    return state_transition_matrix
end
function transitionMatrix( states::Array{T1,2}, state_matrix::Array{T2,2}, nb_type::T4, type_series::Vector{T5}, min_lag::T6, max_lag::T7, d_lag::T8) where { T1 <: Int, T2 <: Real, T4 <: Real, T5 <: Real, T6 <: Real, T7 <: Int, T8 <: Int }

    nb_states=size(states)[1]
    nb_series=size(state_matrix)[1]
    nb_steps=size(state_matrix)[2]
    nb_lag_points=Int(trunc((max_lag-min_lag)/d_lag))

    state_transition_matrix=zeros( nb_states, nb_states, nb_lag_points )
    count_lag=1
    for lag=min_lag:d_lag:max_lag-1
        print("Computing Transition Matrix - Progress: ",lag/max_lag*100,"%\n")
        for serie=1:nb_series
            for j=lag+1:nb_steps
                if state_matrix[serie,j-lag] == 0 || state_matrix[serie,j] == 0
                    continue
                end
                state_transition_matrix[ state_matrix[serie,j-lag], state_matrix[serie,j], count_lag ] += 1
            end
        end
        count_lag += 1
    end

    # Normalization
    for lag=1:nb_lag_points
        print("Normalizing Transition Matrix: ",lag/nb_lag_points*100,"%\n")
        for i=1:nb_states
            sum_transition=sum( state_transition_matrix[i,:,lag] )
            if sum_transition != 0
                state_transition_matrix[i,:,lag] /= sum_transition
            end
        end
    end

    return state_transition_matrix
end
function chappmanKormologov( transition_matrix::Array{T1,3} ) where { T1 <: Real }
    nb_states   = size(transition_matrix)[1]
    nb_lag_time = size(transition_matrix)[3]
    nb_lag_compare = Int(trunc(nb_lag_time))
    transition_matrix_kolmo= zeros(nb_states,nb_states,nb_lag_compare)
    for lag=1:nb_lag_compare
        print("Chappman Kolmogorov Test - ",lag/nb_lag_compare*100,"%\n")
        for i=1:nb_states
            for j=1:nb_states
                for k=1:nb_states
                    transition_matrix_kolmo[i,j,lag] += transition_matrix[i,k,lag]*transition_matrix[k,j,lag]
                end
            end
        end
    end
    return transition_matrix_kolmo
end
function writeStates( file::T1 , states::Array{T2,2}, count_states::Vector{T3}, types::Vector{T4} ) where { T1 <: AbstractString, T2 <: Real, T3 <: Real , T4 <: AbstractString }
    file_out=open(file,"w")
    n_dim = size( states)[2]
    nb_states=size(states)[1]
    nb_types=size(types)[1]
    for i=1:nb_types
        for j=1:Int(n_dim/nb_types)
            write(file_out,string(types[i]," "))
        end
    end
    write(file_out,string("\n"))
    for i=1:nb_states
        for j=1:n_dim
            write(file_out,string(Int(states[i,j])," "))
        end
        write(file_out,string(count_states[i],"\n"))
    end
    close(file_out)
    return
end
function writeStateMatrix( file::T1, state_matrix::Array{T2,2}) where { T1 <: AbstractString, T2 <: Real }
    nb_steps=size(state_matrix)[2]
    nb_series=size(state_matrix)[1]
    file_out=open(file,"w")
    for step=1:nb_steps
        for serie=1:nb_series
            write(file_out,string(state_matrix[serie,step]," "))
        end
        write(file_out,string("\n"))
    end
    close(file_out)
    return
end
function writeTransitionsMatrix( file::T1, transitions_matrix::Array{T2,3} ) where { T1 <: AbstractString, T2<: Real}
    nb_states=size(transitions_matrix)[1]
    nb_steps=size(transitions_matrix)[3]

    file_out=open(file,"w")
    for step=1:nb_steps
        write(file_out,string(step," "))
        for state_i=1:nb_states
            for state_f=1:nb_states
                write(file_out,string(transitions_matrix[state_i,state_f,step]," "))
            end
        end
        write(file_out,string("\n"))
    end
    close(file_out)

    return
end
function writeTransitionsMatrix( file::T1, transitions_matrix::Array{T2,3}, file2::T3, transition_matrix_CK::Array{T3,3} ) where { T1 <: AbstractString, T2 <: Real, T3 <: Real }
    nb_states=size(transitions_matrix)[1]
    nb_steps=size(transitions_matrix)[3]

    file_out=open(file,"w")
    for step=1:nb_steps
        write(file_out,string(step," "))
        for state_i=1:nb_states
            for state_f=1:nb_states
                write(file_out,string(transitions_matrix[state_i,state_f,step]," "))
            end
        end
        write(file_out,string("\n"))
    end
    close(file_out)

    nb_states=size(transitions_matrix_CK)[1]
    nb_steps=size(transitions_matrix_CK)[3]

    file_out=open(file2,"w")
    for step=1:nb_steps
        write(file_out,string(step*2," "))
        for state_i=1:nb_states
            for state_f=1:nb_states
                write(file_out,string(transitions_matrix_CK[state_i,state_f,step]," "))
            end
        end
        write(file_out,string("\n"))
    end
    close(file_out)

    return
end
function readStateMatrix( file::T1 ) where { T1 <: AbstractString }
    file_in=open(file)
    lines=readlines(file_in)
    close(file_in)

    nb_steps=size(lines)[1]
    nb_states=size(split(lines[1]))[1]

    state_matrix=zeros(nb_steps,nb_states)
    for step=1:nb_steps
        print("Progress: ",step/nb_steps*100,"%\n")
        line=split(lines[step])
        for states=1:nb_states
            state_matrix[step,states] = parse(Float64,line[states])
        end
    end

    return state_matrix
end
function readStates( file::T1 ) where { T1 <: AbstractString }
    file_in=open(file)
    lines=readlines(file_in)
    close(file_in)

    nb_states=size(lines)[1]-1
    nb_dim=size(split(lines[1]))[1]

    states=zeros(nb_states,nb_dim)
    for state=1:nb_states
        line=split(lines[state+1])
        for i=1:nb_dim
            states[state,i] = parse(Float64,line[i+1])
        end
    end

    return states, nb_states
end
function readTransitionMatrix( file::T1 ) where { T1 <: AbstractString }
    file_in=open(file)
    lines=readlines(file_in)
    close(file_in)

    nb_states=Int(sqrt(size(split(lines[1]) )[1]-1))
    nb_steps=size(lines)[1]

    transition_matrix=zeros(nb_states,nb_states,nb_steps)
    for step=1:nb_steps
        line=split(lines[step])
        for state_i=1:nb_states
            for state_j=1:nb_states
                transition_matrix[state_i,state_j,step] = parse(Float64,line[ (state_i-1)*nb_states+state_j+1 ])
            end
        end
    end

    return transition_matrix
end

end
