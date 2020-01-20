module press_stress

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

function diagStressTensor( stress_tensor_matrix::Array{T1,3} ) where { T1 <: Real }
    nb_step=size(stress_tensor_matrix)[1]
    vals=zeros(nb_step,3)
    for step=1:nb_step
        vals[i,:] = diagStressTensor(stress_tensor_matrix[i,:,:])
    end
    return vals
end

end
