module stats

# Loading necessary modules from LibAtomicSim
using Statistics
using Bootstrap

# Description:
# Set of functions used for stastitics:
# -> Block average
# -> Bootstrap
# -> Histogram
# -> 2D Histogram

# Block averages
#------------------------------------------------------------------------------
# Compute block averages with a given size of blocks
function blockAverage( data::Vector{T1}, size_block::T2 ) where { T1 <: Real, T2 <: Int }
    # Argument
    # - data: Vector of real containing data to block average
    # - size_block: int, size of the blocks
    # Output
    # - Average of the blocks
    # - Variance over the blocks

    # Get the number of data point
    nb_point = size(data)[1]

    # Number of blocks
    nb_block = round(Int, nb_point/size_block )

    # Initialize averages for all blocks
    averages_block=zeros(nb_block)

    # Loop over blocks
    for block=1:nb_block-1
        # Computes averages of data for each block
        averages_block[block] = Statistics.mean( data[ (block-1)*size_block+1 : block*size_block ] )
    end

    # Returns average of the averages and standard deviations over the averages of the blocks
    return Statistics.mean(averages_block), Statistics.var(observation_block)/(n_block-1)
end
# Compute block averages with different sizes for blocks - to check convergence
function blockAverage( data::Vector{T1}, min_block_size::T2, max_block_size::T3, stride_block::T4 ) where { T1 <: Real, T2 <: Int, T3 <: Int, T4 <: Int }
    # Argument
    # - data: Vector of data to block average
    # - min_block_size, max_block_size: Min and max size of the blocks
    # - stride_block: stride to use to go from min to max of block sizes
    # Output
    # - blocks: size of the blocks
    # - meanBlocks: associated block average with the size of blocks
    # - varBlocks: associated block variance with size of blocks

    # Number of data points
    nb_point = size( data )[1]

    # Number of block sizes
    nb_blocks = Int( ( max_block_size-min_block_size )/stride_block ) + 1

    # Initialize outputs
    meanBlocks = zeros( block_nb ) # - Means
    varBlocks  = zeros( block_nb ) # - Variances
    blocks     = zeros( block_nb ) # - Blocks sizes

    # Initialize block size counter at 1
    block_ctrl=1

    # Loop over blocks sizes
    for block_size=min_block_size:stride_block:max_block_size
        # Compute block average for each block size
        meanBlocks[block_ctrl], varBlocks[block_ctrl] = blockAverage( data, block_size )

        # Get the block size
        blocks[block_ctrl] = block_size

        # Number of blocks size looped
        block_ctrl += 1
    end

    # Return size of the blocks, block averages values - average and variance
    return blocks, meanBlocks, varBlocks
end
#-------------------------------------------------------------------------------

# Bootstrap functions
#-------------------------------------------------------------------------------
# Simple bootstrap with external functions
function bootstrap( data::Vector{T1}, n_boot::T2) where { T1 <: Real, T2 <: Int }
    # Argument:
    # - data: vector of data to bootstrap
    # - n_boot: number of replicas
    # Output
    # - Average computed by Bootstrap
    # - Standard deviation computed by Boostrap

    # Compute Bootstrap with Balanced sampling using Bootstrap module
    bs = bootstrap( Statistics.mean, data, Bootstrap.BalancedSampling( n_boot ) )

    # Return average and standard deviation
    return bs.t0[1], stderror( bs )[1]
end
#-------------------------------------------------------------------------------


# Histogram
#-------------------------------------------------------------------------------
# Creates a simple histogram by giving data as vector, the number of box, min and max of the boxes
function histogram( data::Vector{T1}, nb_box::T2, min_hist::T3, max_hist::T4 ) where { T1 <: Real , T2 <: Int, T3 <: Real, T4 <: Real }
    # Argument
    # - data: vector of real, contains the data
    # - nb_box: number of histogram boxes
    # - min_ : minimum of the boxes
    # - max_ : maximum of the boxes
    # Output
    # - histogram: Vector with Int, containing number of occurences of data within each boxes

    # Initialize boxes
    hist = zeros( nb_box )

    # Computes size of the boxes
    delta_box = ( max_hist - min_hist )/nb_box

    # Loop over boxes
    for box=1:size(data)[1]
        # Box_integer:
        box_int = round(Int, ( data[box] - min_hist )/delta_box ) + 1
        # We keep only elements that are valid
        if box_int > 0 && box_int <= nb_box
            hist[ box_int ]  = hist[ box_int ] +  1
        end
        # Adds point to box
    end

    # Return histogram
    return hist
end
# Creates a simple histogram by giving data as vector, the number of box
# - uses min and max of the data
function histogram( data::Vector{T1}, nb_box::T2) where { T1 <: Real , T2 <: Int }
    # Argument
    # - data: vector of real, contains the data
    # - nb_box: number of boxes of the histogram
    # Output
    # - histogram of the data

    # Get min and max of the data
    min_ = minimum(data)
    max_ = maximum(data)

    # Returns the histogram
    return histogram( data, nb_box, min_, max_ )
end
# Creates a simple normed histogram by giving data as vector, the number of box
# - uses min and max of the data
function histogramNormed( data::Vector{T1}, nb_box::T2) where { T1 <: Real , T2 <: Int }
    # Argument
    # - data: Vector of real, containing data
    # - nb_box: Int, number of boxes of the histogram
    # Output
    # - histogram: histogram of the data, normed so that sum(histogram=1)

    # Computes the histogram of the data using min and max of the data
    hist = histogram( data, nb_box )

    # Normalize histogram
    hist = hist./sum(hist)

    # Returns histogram
    return hist
end
# Creates a simple normed histogram by giving data as vector, the number of box
# - uses min and max of the data
function histogramNormed( data::Vector{T1}, nb_box::T2, min_hist::T3, max_hist::T4 ) where { T1 <: Real , T2 <: Int, T3 <: Real, T4 <: Real }
    # Argument
    # - data: Vector of real, containing data
    # - nb_box: Int, number of boxes of the histogram
    # - min_, max_: min and max of the histogram boxes
    # Output
    # - hist: histogram of the data, normed so that sum(histogram=1)

    # Computes the histogram of the data
    hist = histogram( data, nb_box, min_hist, max_hist )

    # Normalize histogram
    hist = hist ./ sum( hist )

    # Returns the histogram
    return hist
end
# Writes histogram to file, using path of the file, assumes uniform histogram
function writeHistogram( file_out::T1, histogram::Vector{T2}, min_::T4, max_::T5 ) where { T1 <: AbstractString, T2 <: Real, T3 <: Int, T4 <: Real, T5 <: Real }
    # Argument
    # - file_out: path of the file to write in
    # - histogram: histogram to write to file
    # - min_, max_ : minimum and maximum center of box the histogram
    # Output
    # - Bool: true if all went well

    # Get the number of boxes of the histogram
    nb_box = size(histogram)[1]

    # Compute the size of the boxes of the histogram
    delta_box = ( max_- min_ )/nb_box

    # Opens file
    handle_out = open( file_out, "w" )

    # Loop over boxes
    for box=1:nb_box
        # Writes data for each box in file
        Base.write( handle_out, string( box*delta_box + min_, " ", histogram[box], "\n" ) )
    end

    # Close file
    close(handle_out)

    # Returns true if all went well
    return true
end
#-------------------------------------------------------------------------------

# 2D histograms
#-------------------------------------------------------------------------------
# Creates a 2D histogram using number of boxes, min and max in each direction
function histogram2D( data::Array{T1,2}, nb_box::Vector{T2}, mins_::Vector{T3}, maxs_::Vector{T4} ) where { T1 <: Real, T2 <: Int, T3 <: Real, T4 <: Real }
    # Arguments
    # - data: 2D real array with data
    # - nb_box: vector of int containing the number of boxes in each dimensions
    # - mins_, maxs_ : vectors of reals containing the min and max in each dimensions
    # Output
    # - 2D histogram: Real 2D array for the histogram

    # Initialize vector for size of boxes
    delta_box = zeros(Int,2)

    # Loop over dimensions
    for dim=1:2
        # Compute the size of the boxes in each dimensions
        delta_box[dim] = round(Int, ( maxs_[dim] - mins_[dim] )/nb_box[dim] )
    end

    # Initialize histogram
    hist = zeros( nb_box[1], nb_box[2] )

    # Loop over data points
    for point=1:nb_point
        # Init box_int vector
        box_int = zeros(Int,2)

        # Determine boxes where to add data point
        # - Loop over dimension
        for i=1:2
            # Compute box for the data
            box_int[i] = round(Int, ( data[point,i] - mins_[i] )/delta_box[i] ) + 1
            # If box number is over the max, put point in the max
            if box_int[i] > 0 && box_int[i] <= nb_box
                hist[ box_int[1], box_int[2] ] = hist[ box_int[1], box_int[2] ] +  1
            end
        end

        # Adds point to the histogram
        hist[ box_int[1], box_int[2] ]  += 1
    end

    # Returns the histogram
    return hist
end
# Creates a 2D histogram using number of boxes
# - uses mins and maxs from the data
function histogram2D( data::Array{T1,2}, nb_box::Vector{T2} ) where { T1 <: Real, T2 <: Int }
    # Arguments
    # - data: 2D array with the data
    # - nb_box: number of boxes in each dimensions
    # Output
    # - A histoggram of the data

    # Initialize vectors for the minimum and maximum
    mins_ = zeros(2)
    maxs_ = zeros(2)

    # Computes the minimum and maximum for each dimension
    for dim=1:2
        mins_[dim] = minimum( data[:,1] )
        maxs_[dim] = minimum( data[:,2] )
    end

    # Returns histogram
    return histogram2D( data, nb_box, mins_, maxs_)
end
# Creates a normed 2D histogram using number of boxes
# - uses mins and maxs from the data
function histogram2DNormed( data::Array{T1,2}, nb_box::Vector{T2} ) where { T1 <: Real, T2 <: Int }
    # Arguments
    # - data: 2D array containing data
    # - nb_box: number of boxes in each dimensions
    # Output
    # - histogram: normed histogram

    # Computes 2D histogram
    hist = histogram2D( data, nb_box )

    # Normalize 2D histogram
    hist /= sum(hist)

    # Returns normalized histogram
    return hist
end
# Creates a normed 2D histogram using number of boxes, mins_ and maxs_ in each direction
function histogram2DNormed( data::Array{T1,2}, nb_box::Vector{T2}, mins_::Vector{T3}, maxs_::Vector{T4} ) where { T1 <: Real, T2 <: Int, T3 <: Real, T4 <: Real }
    # Argument
    # - data : 2D array containing all data to be used
    # - nb_box: number of box in each dimensions
    # - mins_, maxs_: minimum and maximum
    # Output
    # histogram: 2D normed histogram

    # Compute 2D histogram
    hist = histogram2D( data, nb_box, mins_, maxs_ )

    # Normalize histogram
    hist /= sum(hist)

    # Returns histogram
    return hist
end
# Writes 2D histogram to file
function writeHistogram2D( file_out::T1, histogram::Vector{T2}, nb_box::T3, mins_::T4, maxs_::T5 ) where { T1 <: AbstractString, T2 <: Real, T3 <: Int, T4 <: Real, T5 <: Real }
    # Arguments
    # - file_out: path to the output file
    # - histogram: 2D histogram to write to file
    # - nb_box: number of boxes in each dimensions (vector of Int)
    # - mins_,maxs_: minimum and maximum of the data in each dimensions (vectors of Real)
    # Output
    # - Bool: returns true if all went well

    # Initialize vector for
    delta_box = zeros(2)

    # Loop over dimensions
    for i=1:2
        # Compute size of the boxes in each dimensions
        delta_box[i] = ( maxs_[i] - mins_[i] )/nb_box[i]
    end

    # Opens file
    handle_out = open( file_out, "w" )

    # Loop over boxes in dimension 1
    for box1=1:nb_box[1]
        # Loop over boxes in dimension 2
        for box2=1:nb_box[2]
            # Writes data to file
            # - Center of box in dimension 1
            Base.write( handle_out, string( box1*delta_box[1] + mins_[1], " ") )
            # - Center of box in dimension 2
            Base.write( handle_out, string( box2*delta_box[2] + mins_[2], " ") )
            # - Value of the histogram
            Base.write( handle_out, string( histogram[box1,box2], "\n") )
        end
    end

    # Close output file
    close(handle_out)

    # Returns true if all went well
    return true
end
#-------------------------------------------------------------------------------

end
