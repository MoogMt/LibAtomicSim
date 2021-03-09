module utils

# Description
#-------------------------------------------------------------
# Various functions that are very useful but that do not really fit very well
# anywhere else
# - basic I/O handling files:
# - basic String handling
# - vector and matrix checking and constructions and swaps
# - switching functions constructions
# - gaussian functions
# - creates point clouds
#-------------------------------------------------------------s

# Export functions
export writeData, writeBasicData
export getNbLines
export getAllFilesWithExtension, getFileName, getFilesName
export copyLine2file, skipLines, readParseLine
export ifLongerCut, char2int, spaces
export nbStepStriding, strideData!
export checkDimVec, checkMatDim, sequenceMatrixH, simpleSequence, swap!
export switchingFunction
export gauss
export createBlobs, createRing

# Basic IO handling
#-------------------------------------------------------------------------------
# Write 1D data to file
function writeData( file_path::T1, data::Vector{T2} ) where { T1 <: AbstractString, T2 <: Real }
    # Argument
    # - file_path: path to the output file
    # - data: vector of data
    # Output
    # - Bool: whether the writting was successful

    # Number of data point
    nb_data=size(data)[1]

    # Open output file
    file_out=open(file_path,"w")

    # Loop over data points
    for point=1:nb_data
        # Writes data to file
        write( file_out, string( data[point], "\n" ) )
    end

    # Closing file
    close(file_out)

    # Sends true if it worked
    return true
end
# Write a vector of x and f(x) data into a file
function writeData( file_out::T1, dx::Vector{T2}, fx::Vector{T3} ) where { T1 <: AbstractString, T2 <: Real, T3 <: Real }
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
# Write a vector of x and 2D data array into a file
function writeData( file_out::T1, data::Array{T2,2} ) where { T1 <: AbstractString, T2 <: Real }
    # Argument
    # - file_out: path to the output file
    # - data: data-array (nb_data,n_dim) with the data
    # Output
    # - True: If the writting went ok

    # Opens the file
    handle_out = open( file_out, "w" )

    # Get number of data
    nb_data = size(data)[1]

    # Get the dimensional of data
    n_dim = size(data)[2]

    # Loop over data points
    for step=1:nb_data
        # Loop over dimensions
        for dim=1:n_dim
            # Writting data to file
            write( handle_out, string( data[step,dim], " " ) )
        end

        # Write end of line to file
        write( handle_out, string("\n") )
    end

    # Close the file
    close(handle_out)

    # Return true if all went well
    return true
end
# Get all lines from a file
function getLines( file_path::T1 ) where { T1 <: AbstractString }
    # Argument:
    # - file_path: path of the input file
    # Output
    # - Vector of AbstractString containing all lines in the file as string

    # Opens file
    file_in = open( file_path )

    # Reads all lines
    lines = readlines( file_in )

    # Close file
    close( file_in )

    # Returns the lines of the file
    return lines
end
# Get the number of lines in a file, returns a user message and false if the file does not exists
function getNbLines( file_path::T1 ) where { T1 <: AbstractString }
    # Argument
    # - file_path: path to the input file
    # Output
    # - Number of lines of files
    # OR false if something went wrong (with a user message)

    # Check that the file exists
    if ! isfile( file_path )
        # If not, sends a message and returns false
        print("No file found at: ",file_path,"\n")
        return false
    end

    # Init line counter at 0
    count_ = 0

    # Opens file
    handle_in = open( file_path )

    # Loop as long as you can read the file
    while ! eof( handle_in )
        # Read one line
        test=readline( handle_in )
        # Add one line to counter
        count_ += 1
    end

    # Close the file
    close( handle_in )

    # Returns the number of lines in the file
    return count_
end
# Returns all files with a given extension in a given folder
function getAllFilesWithExtension( folder_path::T1, extension::T2 ) where { T1 <: AbstractString, T2 <: AbstractString }
    # Argument
    # - folder_path: path to the target folder
    # - extension: type of the extension to select
    # Output
    # - files: Vector of String with names of files with target extension

    # Check that the directory exists
    if ! isdir( folder_path )
        # If the folder does not exists, return false
        print("Folder ",folder_path," does not exists!\n")
        return false
    end

    # Reads the folder for all the files
    all_files = readdir( folder_path )

    # Initialize vector string
    files = Vector{ AbstractString }( undef, 0 )

    # Loop over files
    for file=1:size(all_files)[1]
        #  Parse file name with "." as deliminator
        keyword = split( all_files[file], "." )

        # Gets number of keys with parsing
        nb_keys = size( keyword )[1]

        # If the last key is the same as the target extension, add the files to the vector
        if keyword[ nb_keys ] == extension
            push!( files, all_files[file] )
        end
    end

    # Return the list of files
    return files
end
# Returns the name of the file, without the extension
function getFileName( file::T1 ) where { T1 <: AbstractString }
    # Argument
    # - file: path to the target file
    # Output
    # - Get the name of the file, minus the extension

    # Parse the file with "." as deliminators
    keyword = split( file, "." )

    # Gets number of keys generated with parsing
    nb_keys = size( keyword )[1]

    # Generate new parser as "." + last key generated
    extension = string( ".", keyword[nb_keys] )

    # Parse with extension as deliminator
    keyword2 = split( file, extension )

    # Return the name of the file
    return keyword2[1]
end
# Get the names of all files, without the extension
function getFilesName( files::Vector{T1} ) where { T1 <: AbstractString }
    # Argument:
    # - files: vector of strings containing the names of the files to parse
    # Output
    # - names: the names of the files to parse

    # Initialize vector of names
    names = Vector{ AbstractString }( undef, size(files)[1] )

    # Loop over the list of file names
    for i=1:size(files)[1]
        # Parse each file
        names[i] = getFileName( files[i] )
    end

    # Return the names
    return names
end
# Reads a given number of lines without getting the data
function skipLines( handle::T1, nb_line::T2 ) where { T1 <: IO, T2 <: Int  }
    # Argument
    # - handle: handler of the input file
    # - nb_line: number of line to skip

    # Loop over the number of line to skip
    for i=1:nb_line
        # Reading line, results goes into emptyness
        test=readline( handle )
    end

    # Returns true if all went well
    return true
end
# Read and parse line with " " deliminator
function readParseLine( file_io::T1 ) where { T1 <: IO }
    # Argument
    # - file_io: input file handler
    # Output
    # - Vector of string with all the parsed element from the line

    # Reads and parse the line with the target deliminator
    return split( readline(file_io) )
end
# Read and parse a line with a given deliminator
function readParseLine( file_io::T1, deliminator::T2 ) where { T1 <: IO, T2 <: AbstractString }
    # Argument
    # - file_io: input file handler
    # - deliminator: string that will be used for parsing
    # Output
    # - Vector of string with all the parsed element from the line

    # Reads and parse the line with the target deliminator
    return split( readline(file_io), deliminator )
end
#-------------------------------------------------------------------------------

# Strings Handling
#-------------------------------------------------------------------------------
# Cuts a string corresponding to a number if it is longer than a given number of characters
function ifLongerCut( string_::T1, length_max::T2, cut_ind::T3 ) where { T1 <: AbstractString, T2 <: Int, T3 <: Int }
    # Argument
    # - string_: string to potentially cut
    # - length_max: maximum length of the string
    # - cut_ind: maximum precision on the number
    # Output
    # - string_ : potentially cut string

    # Check if the string is longer than a given amount
    if length(string_) > length_max
        # parse the string into a number
        number = round( parse(Float64,string_), digits=cut_ind )

        # Translates number back into a string
        string_ = string_(number)
    end

    # Returns the potentially cut string
    return string_
end
# Transforms a character into a int
function char2int(char::T) where { T <: Char }
    # Argument
    # - char: character to turn into int
    # Output
    # - Int that corresponds to the character

    # All characters in CapsLock
    maj=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
    # All characters in normal size
    min=['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']

    # Loop over characters, do not respect casse
    for i=1:26
        # If the character in the array matches the target one, we replace it by the associated index
        if char == maj[i] || char == min[i]
            return i
        end
    end
end
# Adds a given number of strings to a string
function spaces( string1::T1, nb::T2 ) where { T1 <: AbstractString, T2 <: Int }
    # Argument
    # - string1: the string to add spaces to
    # - nb: number of spaces to add
    # Output
    # - string2: the string with added spaces

    # Loop over the number of spaces to add
    for i=1:nb
        # Add space to string
        string1 = string( string1, " " )
    end

    # returns the modified string
    return string1
end
#-------------------------------------------------------------------------------

# Striding stuff
#-------------------------------------------------------------------------------
# Returns the number of steps that would results of a given striding
function nbStepStriding( nb_step::T1 , stride_::T2 ) where { T1 <: Int, T2 <: Int }
    # Argument
    # - nb_step: number of steps
    # - stride_: size of the stride
    # Output
    # - Int, number of step remaining if stride_ is used

    # Returns the number of strides is stride_ is module of nb_step
    if nb_step % stride_ == 0
        return Int( nb_step/stride_ )
    # Otherwise it's (nb_step/stride_)+1
    else
        return trunc( Int, nb_step/stride_ ) + 1
    end
end
# Returns the number of steps that would results of a given striding, if we ignore the first nb_ignored steps
function nbStepStriding( nb_step::T1, stride_::T2, nb_ignored::T3 ) where { T1 <: Int, T2 <: Int, T3 <: Int }
    # Argument
    # - nb_step: number of steps
    # - stride_: size of the stride
    # - nb_ignored : number of steps to ignore
    # Output
    # - Int, number of step remaining if stride_ is used with nb_ignored steps skipped

    # Returns the number of remaining steps
    return nbStepStriding( nb_step - nb_ignored, )
end
# Striding data
function strideData!( data::Vector{T1}, stride::T2 ) where { T1 <: Real, T2 <: Int }
    # Argument
    # - data: vector with data to stride
    # Output
    # - Vector with a stride
    # Or False if the stride is negative or too large

    # If stride is problematic return false
    if stride < 0 || stride > size(data)[0]
        return false
    # Or return the data with stride
    else
        return data[1:stride:size(data)[0]]
    end
end
#-------------------------------------------------------------------------------

# Vectors
#-------------------------------------------------------------------------------
# Checks whether vector is of dimension dim
function checkDimVec( vector::Vector{T1}, dim::T2 ) where { T1 <: Real, T2 <: Int }
    # Argument
    # - vector: Vector to check the dimension
    # - dim: expected dimension of the vector
    # Output
    # Bool: whether or not the vector has the expected dimension

    # Size of the vector
    sizevec = size(vector)[1]

    # If the vector is the write size, return true
    if ( sizevec == dim )
        return true
    # If the vector is not the wrong size, return false and sends message
    else
        error("Error! Wrong dimension for vector ($sizevec instead of $dim)")
        return false
    end
end
# Checks whether matrix is of dimension (xdim,ydim)
function checkMatDim( matrix::Matrix{T1}, xdim::T2, ydim::T3 ) where { T1 <: Real, T2 <: Int, T3 <: Int }
    # Argument
    # - matrix: 2D matrix of dimension (xdim,ydim)
    # - xdim, ydim: matrix dimensions
    # Output
    # - Bool: whether the matrix has the expected dimensions

    # Actual dimensions of the matrix
    sizemat_x = size(matrix)[1] # X dim
    sizemat_y = size(matrix)[2] # Y dim

    # Check if the matrix has the expected dimension
    if ( sizemat_x == xdim && sizemat_y == ydim )
        # If so, returns true
        return true
    else
        # If not, returns an error message and returns false
        error("Error! Wrong dimension for matrix!\n Got ($sizexmat,$sizeymat) instead of ($xdim,$ydim)")
        return false
    end
end
# Creates a matrix where element follow a simple sequence
function sequenceMatrixH( nb_element::T1 ) where { T1 <: Int }
    # Argument
    # - nb_element: maximum size of the cell (and of the sequence)
    # Output
    # - matrix: (Int,nb_element,nb_element)

    # Initialize matrix with 0
    matrix = zeros( Int, nb_element, nb_element )

    # Loop over dimension x
    for i=1:nb_element
        # Loop over dimension y
        for j=1:nb_element
            # Affect matrix element
            matrix[i,j] = j
        end
    end

    # Returns the sequence matrix
    return matrix
end
# Simple sequence vector creation
function simpleSequence( size_::T1 ) where { T1 <: Int }
    # Argument
    # - size_ : size of the sequence
    # Output:
    # - vector: vector with the simple Sequence

    # Initialize vector
    vector = zeros(Int, size_ )

    # Loop over the size of the sequence
    for point=1:size_
        # Each point of the vector has the index value
        vector[point] = point
    end

    # Return vector
    return vector
end
# Swap in a vector
function swap!( vector::Vector{T1}, index1::T2, index2::T3 ) where { T1 <: Real, T2 <: Int, T3 <: Int }
    # Argument
    # - data : vector (nb_point)
    # - index1, index2: index of data to swap
    # Output
    # - None

    # Storing data
    stock = vector[index1]

    # Moves data_index2 into data_index1
    vector[index1] = vector[index2]

    # Moves data from storage into data_index2
    vector[index2] = stock

    # Returns nothing
    return
end
# Swap in an array
function swap!( array::Array{T1,2}, index1::T2, index2::T3 ) where { T1 <: Real, T2 <: Int, T3 <: Int }
    # Argument
    # - data : vector (nb_point)
    # - index1, index2: index of data to swap
    # Output
    # - None

    # Storing data
    stock = array[index1,:]

    # Moves data_index2 into data_index1
    array[index1,:] = array[index2,:]

    # Moves data from storage into data_index2
    array[index2,:] = stock

    # Returns nothing
    return
end
#-------------------------------------------------------------------------------

# Switching Functions
#-------------------------------------------------------------------------------
# Returns the result of a switching function
function switchingFunction( x::T1, d::T2, n::T3, m::T4) where { T1 <: Real, T2 <: Real, T3 <: Int, T4 <: Int }
    # Argument
    # - x: value where to evaluate the SF
    # - d, n, m: parameters of the switching function
    # Output
    # Results of the switchin function

    # Returns the result of the switching function
    return ( 1 - (x/d)^n )/( 1 - (x/d)^m )
end
# Returns the results of a switching function, in case where m=2*n
function switchingFunction( x::T1, d::T2, n::T3 ) where { T1 <: Real, T2 <: Real, T3 <: Int }
    # Argument
    # - x : point where we want to evaluate the switching function
    # - d, n: parameters of the switching function
    # Output
    # Results of the switching function (scalar)

    # Returns the value
    return 1/( 1 + (x/d)^n )
end
#-------------------------------------------------------------------------------

# Gaussian related functions
#-------------------------------------------------------------------------------
# Returns a gaussian with definite parameters for a single point, with same amplitude, central position and width in all dimensions
function gauss( amplitude::T1, position::Vector{T2}, width::T3,  x :: Vector{T4} ) where { T1 <: Real, T2 <: Real, T3 <: Real, T4 <: Real }
    # Argument
    # - amplitude: amplitude of the gaussian
    # - position: vector with the central position of the gaussian
    # - width: width of the gaussian (same in all directions)
    # - x : position where we want to evaluate the gaussian
    # Output
    # - value of the gaussian at a given point

    # Initialize value
    value=0

    # Loop over vector dimension
    for i=1:size(position)[1]
        # Adds square values around central position
        dist = x[i] - position[i]
        value += dist*dist
    end

    # Return the value of the gaussian
    return amplitude*exp( - ( value )/(2*( width*width ) ) )
end
# Returns a gaussian with definite parameters for a single point, with different amplitude, central positions and widths in different dimensions
function gauss( amplitudes::Vector{T1}, positions::Array{T2,2}, widths::Vector{T3},  x :: Vector{T4} ) where { T1 <: Real, T2 <: Real, T3 <: Real, T4 <: Real }
    # Argument
    # - amplitude: vector of reals with the amplitudes of the gaussian
    # - position: vector with the central position of the gaussian
    # - widths: vector of real for the width of the gaussian (same in all directions)
    # - x : position where we want to evaluate the gaussian
    # Output
    # - value of the gaussian at a given point

    # Initialize value
    value=0

    # Loop over dimensions
    for i=1:size(amplitudes)[1]
        # Add the value of the gaussian in all dimension
        value += gauss( amplitudes[i], positions[i,:], widths[i], x )
    end

    # Returns the value
    return value
end
# Returns a gaussian with definite parameters for several points, with different widths in different dimensions
function gauss( amplitudes::Vector{T1}, positions::Array{T2,2}, widths::Vector{T3},  x :: Array{T4,2} ) where { T1 <: Real, T2 <: Real, T3 <: Real, T4 <: Real }
    # Argument
    # - amplitude: vector of reals with the amplitudes of the gaussian
    # - position: vector with the central position of the gaussian
    # - widths: vector of real for the width of the gaussian (same in all directions)
    # - x : position where we want to evaluate the gaussian
    # Output
    # - values of the gaussian at the target points

    # Number of data points
    nb_points=size(x)[1]

    # Initialize vector for results
    values=zeros(nb_points)

    # Loop over points
    for i=1:nb_points
        # Compute gaussian values for each point
        values[i] = gauss( amplitudes, positions, widths,x[i,:] )
    end

    # Returns vector with values
    return values
end
#-------------------------------------------------------------------------------

# Point cloud shape generation
#-------------------------------------------------------------------------------
# Generates points in spherical blobs
function createBlobs( n_points::Vector{T1} , centers::Array{T2,2}, spread::Vector{T3} ) where { T1 <: Int, T2 <: Real, T3 <: Real }
    # Argument
    # - n_points: number of points in each blobs (vector, int)
    # - centers: centers of the blobs      ( real, array ( nb_blobs, nb_dim ) )
    # - spread: spread of the various blob ( real, array ( nb_blobs, nb_dim ) )
    # Output:
    # - points: positions of the points generated in all blobs

    # Get number of blocs
    n_blobs = size( n_points )[1]

    # Number of total points
    n_points_total = sum( n_points )

    # Dimension of the space
    n_dim = size(centers)[2]

    # Vector with points
	points=zeros( n_points_total, n_dim )

    # Compute square of the spreads of the blobs
	spread_2 = spread.*spread

    # Loop over blobs
	for blob=1:n_blobs
        # Compute the offset for the index of the position of the points
		start_count=sum(n_points[1:blob-1])

        # Loop over the number of points of the blob
		for point=1:n_points[i]
            # First try to put the point somewhere
			try_ = ( rand( n_dim ) .- 0.5 ) * 2 * spread[blob]

            # Loop as long as the point is not within the spread of the blobs
			while sum( try_ .* try_ ) > spread_2[blob]
				try_ = ( rand( n_dim ) .- 0.5 ) * 2 * spread[blob]
			end

            # Once point is successfully positioned, adds position to the list
			points[ start_count + point, :] = centers[blob,:] .+ try_
		end
	end

    # Returns the array of all the points of the blobs
	return points
end
# Generates points in a ring shape
function createRing( n_points::T1, centers::Vector{T2}, small_radius::T3, width::T4 ) where { T1 <: Int, T2 <: Real, T3 <: Real, T4 <: Real }
	n_dim=size(centers)[1]
	points=zeros(Real, n_points,n_dim)
	R=0
	angles=ones(Real, n_dim-1)
	for i=1:n_points
		# Randomize a distnace to the center
		R = small_radius+rand()*width
		angles=rand(n_dim-1)*pi
		angles[n_dim-1]=rand()*2*pi
		for j=1:n_dim-1
			points[i,j]=R*cos(angles[j])
		end
		points[i,n_dim]=R
		for j=2:n_dim
			for k=1:j-1
				points[i,j] = points[i,j]*sin(angles[k])
			end
		end
	end
	return points
end
#-------------------------------------------------------------------------------

end
