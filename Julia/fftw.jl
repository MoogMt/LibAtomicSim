module fftw

using FFTW

function doFourierTransform( signal::Vector{T1}, dt::Real ) where { T1 <: Real, T2 <: Real }
    freq=fftfreq(size(signal)[1],1/dt)
    intensity= abs.(fft(signal))
    return  freq[ freq .> 0 ], intensity[ freq .> 0 ]
end

function doFourierTransformShift( signal::Vector{T1}, dt::T2 ) where { T1 <: Real, T2 <: Real }
    freq=fftfreq(size(signal)[1],1/dt) |> fftshift
    intensity= abs.(fft(signal)) |> fftshift
    return  freq[ freq .> 0 ], intensity[ freq .> 0 ]
end

function doDCT( signal::Vector{T1}, dt::T2 ) where { T1 <: Real, T2 <: Real }
    freq=fftfreq(size(signal)[1],1/dt)
    intensity=abs.(dct(signal))
    return  freq[ freq .> 0 ], intensity[ freq .> 0 ]
end

function doDCTShift( signal::Vector{T1}, dt::T2 ) where { T1 <: Real, T2 <: Real }
    freq=fftfreq(size(signal)[1],1/dt) |> fftshift
    intensity=abs.(dct(signal)) |> fftshift
    return  freq[ freq .> 0 ], intensity[ freq .> 0 ]
end

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
