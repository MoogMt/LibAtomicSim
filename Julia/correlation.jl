module correlation

# Loading necessary modules from default repository
using Statistics

# Description
# - Sets of functions dealing with correlation/autocorrelation of
# a given signal

# Export functions
export autocorr, autocorrAvg, autocorrAvgSig, autocorrNorm

# Autocorrelation functions
#-------------------------------------------------------------------------------
function autocorr( signal::Vector{T1}, max_lag::T2 ) where { T1 <: Real, T2 <: Int }
    # Argument
    # - signal: signal to compute the aucorrelation of (vector real)
    # - max_lag: maximum lag of the autocorrelation
    # Output
    # - autocor_signal: Autocorrelation of the signal

    # Get number of steps of signal
    nb_step=size(signal)[1]

    # Initialize vector for autocorrelation of signal
    autocor_signal=zeros(max_lag)

    # Loop over tau
    for lag=0:max_lag-1
        # Loop over steps
        for step=1:nb_step-lag
            # Computes autocorrelation for a given tau
            autocor_signal[ lag + 1 ] += ( signal[step] )*( signal[ step + lag ] )
        end

        # Normalize autocorrelation step
        autocor_signal[ lag + 1 ] /= ( nb_step - lag )
    end

    # Returns autocorrelation signal
    return autocor_signal
end
# Autocorrelation of the signal normalized by its average
function autocorrAvg( signal::Vector{T1}, max_lag::T2 ) where { T1 <: Real, T2 <: Int }
    # Argument
    # - signal: signal to compute the aucorrelation of (vector real)
    # - max_lag: maximum lag of the autocorrelation
    # Output
    # - autocor_signal: Autocorrelation of the signal

    # Get number of steps of signal
    nb_step = size(signal)[1]

    # Initialize vector for autocorrelation of signal
    autocor_signal = zeros(Real, max_lag )

    # Computes average of the signal
    avg_sig = mean( signal )

    # Loop over tau
    for lag=0:max_lag-1
        # Loop over steps
        for step=1:nb_step-lag
            # Computes autocorrelation
            autocor_signal[lag+1] += ( signal[step] - avg_sig )*( signal[step+lag] - avg_sig )
        end

        # Computes normalization factor
        autocor_signal[lag+1] /= (nb_step-lag)
    end

    # Returns autocorrelation signal
    return autocor_signal
end
# Autocorrelation of the signal normalized by its average and variance
function autocorrAvgSig( signal::Vector{T1}, max_lag::T2 ) where { T1 <: Real, T2 <: Int }
    # Argument
    # - signal: signal to compute the aucorrelation of (vector real)
    # - max_lag: maximum lag of the autocorrelation
    # Output
    # - autocor_signal: Autocorrelation of the signal

    # Get number of steps of signal
    nb_step=size(signal)[1]

    # Initialize vector for autocorrelation of signal
    autocor_signal=zeros(max_lag)

    # Computes average of the signal
    avg_sig=mean(signal)

    # Computes average of the signal
    var_sig=Statistics.var(signal)

    # Loop over tau
    for lag=0:max_lag-1
        # Loop over steps
        for step=1:nb_step-lag
            # Computes autocorrelation of the signal for a given lag
            autocor_signal[lag+1] += ( signal[step] - avg_sig)*( signal[step+lag] - avg_sig )
        end

        # Normalize autocorrelation signal
        autocor_signal[lag+1] /= ( ( nb_step - lag )*var_sig )
    end

    # Returns autocorrelation of signal
    return autocor_signal
end
# Computes autocorrelation of signal, and then normalize with regard to the first step
function autocorrNorm( signal::Vector{T1}, max_lag::T2 ) where { T1 <: Real, T2 <: Int }
    # Argument
    # - signal: vector with the signal (vector real)
    # - max_lag: maximum lag of the signal
    # Output
    # - aucocor_signal: autocorrelation of the signal

    # Computes autocorrelation of the signal
    autocor_signal = autocorr( signal, max_lag )

    # Normalize by the first element of the vector
    autocor_signal /= autocor_signal[1]

    # Returns the normalized autocorrelation signal
    return autocor_signal
end
#-------------------------------------------------------------------------------

# EXEMPLES OF USE
#-------------------------------------------------------------------------------
# x=range(0,stop=2*pi,step=0.005)
# y=sin.(x)
# z=correlation.autocorrNorm(y,Int(trunc(size(y)[1]*0.8)))
#
# file_out=open(string("/home/moogmt/test.dat"),"w")
# for i=1:size(z)[1]
#     Base.write(file_out,string(i," ",z[i],"\n"))
# end
# close(file_out)
#-------------------------------------------------------------------------------

end
