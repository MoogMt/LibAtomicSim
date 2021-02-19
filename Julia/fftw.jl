module fftw

# Loading necessary module from Julia default repository
using FFTW

# Description
# - Allows easier computations of FFTW and DCT

# Export all functions
export doFourierTransform, doFourierTransformShift
export doDCT, doDCTShift

# Computes Fast Fourrier Transform
#-------------------------------------------------------------------------------
# Simple FFT + Associated functions
function doFourierTransform( signal::Vector{T1}, dt::Real ) where { T1 <: Real, T2 <: Real }
    # Argument
    # - signal: values of the signal (real, vector, nb_step)
    # - dt: timestep of the signal
    # Output
    # - freq: vector of real with frequencies associated with FFT
    # - intensity: intensity of the FFT

    # Computes frequencies
    freq = fftfreq( size( signal )[1], 1/dt )

    # Computes FFTW of the signal
    intensity = abs.( fft( signal ) )

    # Returns frequencies of FFT and FFT of the signal
    return  freq[ freq .> 0 ], intensity[ freq .> 0 ]
end
function doFourierTransformShift( signal::Vector{T1}, dt::T2 ) where { T1 <: Real, T2 <: Real }
    # Argument
    # - signal: values of the signal (real, vector, nb_step)
    # - dt: timestep of the signal
    # Output
    # - freq: vector of real with frequencies associated with FFT
    # - intensity: intensity of the FFT

    # Computes frequencies
    freq = fftfreq( size( signal )[1], 1/dt ) |> fftshift

    # Computes FFTW of the signal
    intensity= abs.(fft(signal)) |> fftshift

    # Returns frequencies of FFT and FFT of the signal
    return  freq[ freq .> 0 ], intensity[ freq .> 0 ]
end
#-------------------------------------------------------------------------------

# Computes Discrete Cosine Transform
#-------------------------------------------------------------------------------
# Simple DCT + associated frequency
function doDCT( signal::Vector{T1}, dt::T2 ) where { T1 <: Real, T2 <: Real }
    # Argument
    # - signal: values of the signal (real, vector, nb_step)
    # - dt: timestep of the signal
    # Output
    # - freq: vector of real with frequencies associated with DCT
    # - intensity: intensity of the DCT

    # Computes frequencies
    freq=fftfreq(size(signal)[1],1/dt)

    # Computes DCT of the signal
    intensity=abs.(dct(signal))

    # Returns frequencies of DCT and DCT of the signal
    return  freq[ freq .> 0 ], intensity[ freq .> 0 ]
end
# Simple DCT + associated frequency shifted to > 0
function doDCTShift( signal::Vector{T1}, dt::T2 ) where { T1 <: Real, T2 <: Real }
    # Argument
    # - signal: values of the signal (real, vector, nb_step)
    # - dt: timestep of the signal
    # Output
    # - freq: vector of real with frequencies associated with DCT
    # - intensity: intensity of the DCT

    # Computes frequencies
    freq      = fftfreq( size( signal )[1], 1/dt ) |> fftshift

    # Computes DCT of the signal
    intensity = abs.( dct( signal ) ) |> fftshift

    # Returns frequencies of DCT and DCT of the signal
    return  freq[ freq .> 0 ], intensity[ freq .> 0 ]
end
#-------------------------------------------------------------------------------

end

#  Examples
# dt=0.001
# x=range(0,stop=5*pi,step=dt)
# y=cos.(2*pi*x*50)
#
# j,k=doFourierTransform(y,dt)
#
# file_out=open(string("/home/moogmt/test2.dat"),"w")
# for i=1:size(z)[1]
#     Base.write(file_out,string(x[i]," ",y[i],"\n"))
# end
# close(file_out)
#
# file_out=open(string("/home/moogmt/test.dat"),"w")
# for i=1:size(k)[1]
#     Base.write(file_out,string(j[i]," ",k[i],"\n"))
# end
# close(file_out)
