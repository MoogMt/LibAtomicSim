module utils

#-------------------------------------------------------------------------------
#  List of useful random functions that I could not fit elsewhere
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

#===============#
# General Array #
#===============#

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

end
