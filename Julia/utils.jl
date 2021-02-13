module utils

# Description
# Various functions that are very useful but that do not really fit very well
# anywhere else
# - basic I/O handling files:
#  -> get number of lines/test existence of file

# TODO: Export all functions

export getNbLines


# Write a vector of x and f(x) data into a file
function writeBasicData( file_out::T1, dx::Vector{T2}, fx::Vector{T3} ) where { T1 <: AbstractString, T2 <: Real, T3 <: Real }
    # Argument
    # - file_out: path to the output file
    # - dx: Input of the function (vector Real; (nb_data))
    # - fx: value associated to dx (vector Real; (nb_data))
    # Output
    # - True: If the writting went ok

    # Opens the file
    handle_out = open( file_out, "w" )

    # Get number of data
    nb_data = size(dx)[1]

    # Loop over data points
    for data=1:nb_data
        # Writting data to file
        write( handle_out, string( dx[data], " ", fx[data], "\n" ) )
    end

    # Close the file
    close(handle_out)

end

function getNbLines( file_path::T1 ) where { T1 <: AbstractString }
    if ! isfile( file_path )
        print("No file found at: ",file_path,"\n")
        return false
    end
    count_ = 0
    handle_in = open( file_path )
    while ! eof( handle_in )
        readline( handle_in )
        count_ += 1
    end
    close(handle_in)
    return count_
end

#-------------------------------------------------------------------------------
function getAllFilesWithExtension( folder_path::T1, extension::T2 ) where { T1 <: AbstractString, T2 <: AbstractString }
    if ! isdir( folder_path )
        print("Folder ",folder_path," does not exists!\n")
        return false
    end
    all_files = readdir( folder_path )
    files = Vector{ AbstractString }( undef, 0 )
    for file=1:size(all_files)[1]
        keyword = split( all_files[file], "." )
        nb_keys = size( keyword )[1]
        if keyword[ nb_keys ] == extension
            push!( files, all_files[file] )
        end
    end
    return files
end
function getFileName( file::T1 ) where { T1 <: AbstractString }
    keyword = split( file, "." )
    nb_keys = size( keyword )[1]
    extension = string( ".", keyword[nb_keys] )
    keyword2 = split( file, extension )
    return keyword2[1]
end
function getFilesName( files::Vector{T1} ) where { T1 <: AbstractString }
    names = Vector{ AbstractString }( undef, size(files)[1] )
    for i=1:size(files)[1]
        names[i] = getFileName( files[i] )
    end
    return names
end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
function ifLongerCut( string_::T1, length_max::T2, cut_ind::T3 ) where { T1 <: AbstractString, T2 <: Int, T3 <: Int }
    if length(string_) > length_max
        number=parse(Float64,string_)
        number=round(number,digits=cut_ind)
        string_=string_(number)
    end
    return string_
end
function skipLines( handle::T1, nb_line::T2 ) where { T1 <: IO, T2 <: Int  }
    for i=1:nb_line
        test=readline( handle )
    end
    return true
end
function getLineElements( file_io::T1 ) where {T1 <: IO }
    return split( readline(file_io) )
end
function copyLine2file( line::T1, file_io::T2 ) where { T1 <: AbstractString, T2 <: IO }
    write(file_io, line )
    write(file_io, string("\n"))
    return true
end
function copyLine2file( line_element::Vector{T1}, file_io::T2 ) where { T1 <: AbstractString, T2 <: IO }
    nb_elements = size(line_element)[1]
    for i=1:nb_elements
        write( file_io, string( line_element[i], " " ) )
    end
    write(file_io, string("\n"))
    return true
end
#-------------------------------------------------------------------------------

# Striding stuff
#-------------------------------------------------------------------------------
function nbStepStriding( nb_step::T1 , stride_::T2 ) where { T1 <: Int, T2 <: Int }
    if nb_step % stride_ == 0
        return Int(nb_step/stride_)
    else
        return trunc(Int,nb_step/stride_)+1
    end
end
function nbStepStriding( nb_step::T1, stride_::T2, nb_ignored::T3 ) where { T1 <: Int, T2 <: Int, T3 <: Int }
    return nbStepStriding( nb_step-nb_ignored, )
end
function strideData!( data::Vector{T1}, stride::T2 ) where { T1 <: Real, T2 <: Int }
    if stride < 0 || stride > size(data)[0]
        return false
        return data[1:stride:size(data)[0]]
    end
end
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
function writeData( file_path::T1, data::Vector{T2} ) where { T1 <: AbstractString, T2 <: Real }
    nb_data=size(data)[1]
    file_out=open(file_path,"w")
    for step=1:nb_data
        write(file_out,string(data[step],"\n"))
    end
    close(file_out)
    return true
end
function getLines( file_path::T1 ) where { T1 <: AbstractString }
    file_in = open( file_path )
    lines = readlines( file_in )
    close( file_in )
    return lines
end
#-------------------------------------------------------------------------------


#==============================================================================#
function determineFolderPath( computers_names::Vector{T1}, paths::Vector{T2} ) where { T1 <: AbstractString, T2 <: AbstractString }
    size_vector = size( computers_names )[1]
    host_name = gethostname()
    for i=1:size_vector
        if computers_names[i] == host_name
            return paths[i]
        end
    end
    return false
end
#==============================================================================#

# VECTORS
#==============================================================================#
# - Checks whether vector is of dimension dim
function checkDimVec(vector::Vector{T1}, dim::T2) where {T1 <: Real, T2 <: Int}
    sizevec = size(vector)[1]
    if ( sizevec == dim )
        return true
    else
        error("Error! Wrong dimension for vector ($sizevec instead of $dim)")
        return false
    end
end
# - Checks whether matrix is of dimension xdim*ydim
function checkMatDim( matrix::Matrix{T1}, xdim::T2, ydim::T3 ) where { T1 <: Real, T2 <: Int, T3 <: Int }
    sizematx=size(matrix)[1]
    sizematy=size(matrix)[2]
    if ( sizematx == xdim && sizematy == ydim )
        return true
    else
        error("Error! Wrong dimension for matrix!\n Got ($sizexmat,$sizeymat) instead of ($xdim,$ydim)")
        return false
    end
end
# - Removes the duplicate elements in a vector
function removeDuplicates( vector::Vector{T1} ) where { T1 <: Real }
    i=1; j=2;
    while i < size(vector)[1]
        while j <= size(vector)[1]
            if vector[i] == vector[j]
                deleteat!(vector,j)
            else
                j=j+1
            end
        end
        i=i+1
        j=i+1
    end
    return vector
end
#==============================================================================#

# STRING
#==============================================================================#
# Parsing
#--------------------------------------------------------------
# Prase a string into a real
function str2rl( string::T ) where { T <: AbstractString }
    return parse(Float64,string)
end
# Parse a string into an Int
function str2int( string::T ) where { T <: AbstractString }
    return parse(Int64,string)
end
#---------------------------------------------------------------
# Conversion
# Transforms a character into a int
function char2int(char::T) where { T <: Char }
    maj=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
    min=['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
    for i=1:26
        if char == maj[i] || char == min[i]
            return i
        end
    end
end
#---------------------------------------------------------------
# String
#----------
# Puts nb spaces at string1 end
function spaces( string1::T1, nb::T2 ) where { T1 <: AbstractString, T2 <: Int }
    string2=string1
    for i=1:nb
        string2=string(string2," ")
    end
    return string2
end
#==============================================================================#


#==============================================================================#
function isIn( element, list )
    for i=1:size(list)[1]
        if element == list[i]
            return true
        end
    end
    return false
end
#==============================================================================#
function sequenceMatrixH( nb_element::T1 ) where { T1 <: Int }
    matrix=zeros(Int,nb_element,nb_element)
    for i=1:nb_element
        for j=1:nb_element
            matrix[i,j]=j
        end
    end
    return matrix
end
#==============================================================================#

# Switching Functions
#==============================================================================#
function switchingFunction( x::T1, d::T2, n::T3, m::T4) where { T1 <: Real, T2 <: Real, T3 <: Int, T4 <: Int}
    return (1-(x/d)^n)/(1-(x/d)^m)
end
function switchingFunction( x::T1, d::T2, n::T3 ) where { T1 <: Real, T2 <: Real, T3 <: Int }
    return 1/(1+(x/d)^n)
end
#==============================================================================#

# GAUSSIAN
#==============================================================================#
function gauss( amplitude::T1, position::Vector{T2}, width::T3,  x :: Vector{T4} ) where { T1 <: Real, T2 <: Real, T3 <: Real, T4 <: Real }
    value=0
    for i=1:size(position)[1]
        value += (x[i]-position[i])*(x[i]-position[i])
    end
    return amplitude*exp( - (value)/(2*(width*width)) )
end
function gauss( amplitudes::Vector{T1}, positions::Array{T2,2}, widths::Vector{T3},  x :: Vector{T4} ) where { T1 <: Real, T2 <: Real, T3 <: Real, T4 <: Real }
    value=0
    for i=1:size(amplitudes)[1]
        value += gauss(amplitudes[i],positions[i,:],widths[i],x)
    end
    return value
end
function gauss( amplitudes::Vector{T1}, positions::Array{T2,2}, widths::Vector{T3},  x :: Array{T4,2} ) where { T1 <: Real, T2 <: Real, T3 <: Real, T4 <: Real }
    nb_points=size(x)[1]
    values=zeros(nb_points)
    for i=1:nb_points
        values[i] = gauss( amplitudes, positions,widths,x[i,:])
    end
    return values
end
#==============================================================================#


function histogram( data::Vector{T1}, min_data::T2, max_data::T3, nb_box_data::T4 ) where { T1 <: Real, T2 <: Real, T3 <: Real, T4 < :Int }
    # Argument
    # - data : Vector of real (nb_data), with the data to put in histogram
    # - min_data: minimum value of the histogram
    # - max_data: maximum value of the histogram
    # - nb_box: number of boxes of the histogram
    # Output
    # - dx: centers of the histogram boxes
    # - hist: normed values of the histogram of the data

    # Number of data points
    nb_data = size(data)[1]

    # Computing the size of the boxes
    delta_data = ( max_data - min_data )/nb_box_data

    # Initialize histogram data
    hist = zeros(Real, nb_box_data + 1 )

    # Initialize boxes center positions
    dx = zeros(Real, nb_box_data + 1 )

    # Computing box centers
    for box=1:nb_box_data
        dx[i] = box*delta_ + min_data
    end

    # Loop over data point
    for i=1:nb_data
        point = round(Int, ( abs(distances_to_plan[i]) - min_ )/delta_ + 1 )
        hist[ point ] = hist[point] + 1
    end

    # Normalization of the histogram
    hist /= sum( hist )

    # Returns the center of boxes and histogram
    return dx, hist
end

end
