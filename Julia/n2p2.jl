module n2p2

export writeN2P2conf

function writeN2P2conf( positions::Array{T1,3}, forces::Array{T2,3}, cell_matrix::Array{T3,2}, atom_types::Vector{T4}, charge::Vector{T5}, charge_system::Vector{T6}, energy::Vector{T7}, file_out::T8 ) where { T1 <: Real, T2 <: Real, T3 <: Real, T4 <: AbstractString, T5 <: Real, T6 <: Real, T7<:Real, T8 <: AbstractString }
    nb_steps,nb_atoms,dim=size(positions)
    file=open(file_out,"w")
    for step=1:nb_steps
        write(file,string("begin\n"))
        write(file,string("comment STEP ",step,"\n"))
        for i=1:3
            write(file,string("lattice "))
            for j=1:3
                write(file,string(cell_matrix[i,j]," "))
            end
            write(file,string("\n"))
        end
        for atom=1:nb_atoms
            write(file,string("atom "))
            for i=1:3
                write(file,string(positions[step,atom,i]," "))
            end
            write(file,string(atom_types[atom]," "))
            write(file,string(charge[atom]," "))
            write(file,string("X "))
            for i=1:3
                write(file,string(forces[step,atom,i]," "))
            end
            write(file,string("\n"))
        end
        write(file,string("energy ",energy[step],"\n"))
        write(file,string("charge ",charge_system[step],"\n"))
        write(file,string("end\n"))
    end
    close(file)
end

end
