module press_stress

using conversion

function computePressure( stress_tensor::Array{T1,2} ) where { T1 <: Real }
    pressure=0
    for i=1:3
        pressure += stress_tensor[i,i]
    end
    return pressure/3
end

function computePressure( stress_tensor_matrix::Array{T1,3} ) where { T1 <: Real }
    nb_step=size(stress_tensor_matrix)[1]
    pressure=zeros(nb_step)
    for step=1:nb_step
        pressure[step] = computePressure( stress_tensor_matrix[step,:,:] )
    end
    return pressure
end

function diagStressTensor( stress_tensor::Array{T1,2} ) where { T1 <: Real }
    return eigvals( stress_tensor )
end

function writePressure( file_path::T1, pressure::Vector{T2} ) where { T1 <: AbstractString, T2 <: Real }
    handle_out = open( file_path, "w" )
    nb_step = size(pressure)[1]
    for step = 1:nb_step
        write( handle_out, string( pressure[step], "\n" ) )
    end
    close( handle_out )
end

function readPressure( file_path::T1 ) where { T1 <: AbstractString }
    if ! isfile( file_path )
        print("No Pressure file at: ",file_path," !\n")
        return false
    end
    handle_in = open( file_path )
    nb_step=0
    while !eof( handle_in )
        readline( handle_in )
        nb_step += 1
    end
    seekstart( handle_in )
    pressure = zeros(nb_step)
    for step=1:nb_step
        pressure[step] = parse( Float64, split( readline( handle_in ) )[1] )*conversion.kbar2Gpa
    end
    close(handle_in)
    return pressure
end

# function diagStressTensor( stress_tensor_matrix::Array{T1,3} ) where { T1 <: Real }
#     nb_step=size(stress_tensor_matrix)[1]
#     vals=zeros(nb_step,3)
#     for step=1:nb_step
#         vals[i,:] = diagStressTensor(stress_tensor_matrix[i,:,:])
#     end
#     return vals
# end

end
