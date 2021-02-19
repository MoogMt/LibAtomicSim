module press_stress

# Loading necessay module from LibAtomicSim
using conversion

# Descriptions
# - functions used to deal with stress tensor and pressure

# Functions
#-------------------------------------------------------------------------------
# Compute pressure from stress tensor
function computePressure( stress_tensor::Array{T1,2} ) where { T1 <: Real }
    # Argument
    # - stress_tensor: matrix containing the stress tensor
    # Output
    # - pressure (real scalar) on the cell

    # Initialize pressure
    pressure=0

    # Compute trace of the stress tensor, loop over dimensions
    for i=1:3
        pressure += stress_tensor[i,i]
    end

    # returns pressure
    return pressure/3
end
# Compute pressure at each step from stress tensor tensor
function computePressure( stress_tensor_matrix::Array{T1,3} ) where { T1 <: Real }
    # Argument
    # - stress_tensor_matrix: tensor containing the stress tensor for each step
    # Output
    # - pressure: vector containing pressure at each time

    # Get the number of steps of the trajectory
    nb_step = size( stress_tensor_matrix )[1]

    # Initialize pressure vector
    pressure=zeros(nb_step)

    # Loop over steps
    for step=1:nb_step
        # Computes pressure at each step
        pressure[step] = computePressure( stress_tensor_matrix[step,:,:] )
    end

    # returns pressure
    return pressure
end#
#-------------------------------------------------------------------------------

# Diagonalize stress tensor
#-------------------------------------------------------------------------------
function diagStressTensor( stress_tensor::Array{T1,2} ) where { T1 <: Real }
    # Argument
    # - stress tensor: matrix containing the stress tensor
    # Output
    # - Eigenvalues of the stress tensor

    # Returns the eigenvalue of the stress tensor
    return eigvals( stress_tensor )
end
#-------------------------------------------------------------------------------

# Write pressure to file
#-------------------------------------------------------------------------------
function writePressure( file_path::T1, pressure::Vector{T2} ) where { T1 <: AbstractString, T2 <: Real }
    # Arguments
    # - file_path: path to the output file
    # - pressure: vector (real) contains the pressure at each step
    # Output
    # - Bool, whether the writting was successful

    # Opens the output file
    handle_out = open( file_path, "w" )

    # Get the number of steps
    nb_step = size(pressure)[1]

    # Loop over steps
    for step = 1:nb_step
        # Write pressure at the current step to file
        write( handle_out, string( pressure[step], "\n" ) )
    end

    # Close output file
    close( handle_out )

    return true
end
#-------------------------------------------------------------------------------

# Read pressure file (basic file)
#-------------------------------------------------------------------------------
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
end#
#-------------------------------------------------------------------------------

end
