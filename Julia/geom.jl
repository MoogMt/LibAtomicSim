module geom

# Loading necessary modules from default Julia repository
using LinearAlgebra

# Description
# Modules that contains useful functions to deal with basic geometric issues

# Export all functions
export distance
export dist2Line, pointFromLine
export angleAlKash

# Compute the distance between two vectors
#-----------------------------------------------------------------------------
function distance( vec1::Vector{T1}, vec2::Vector{T2} ) where { T1 <: Real , T2 <: Real }
    # Argument
    # - vec1, vec2: the two vectors we want to compute the distance from
    # Output
    # - The distance between the two vectors (real, scalar)

    # return distance between the two vectors
    return LinearAlgebra.norm( vec1 - vec2 )
end
#-----------------------------------------------------------------------------


# Distance of a point to a line
#-----------------------------------------------------------------------------
# Computes the distance of a point to a line
function dist2Line( vector_2_line::Vector{T1}, line_vector::Vector{T2} ) where { T1 <: Real, T2 <: Real }
    # Argument:
    # - vector_2_line: vector from the point to the line (vector real)
    # - line_vector: director vector of the line (vector real)

    # Returns the distance of the point to the line
    return abs( dot( vector_2_line, line_vector ) )/LinearAlgebra.norm( line_vector )
end
# Computes the distance of a specific point ot a line
function pointFromLine( point::Vector{T1}, line_vector::Vector{T2} , point_in_line::Vector{T3} ) where { T1 <: Real, T2 <: Real, T3 <: Real }
    # Argument
    # - point: point we want the distance from the line
    # - line_vector: vector director of the line
    # - point_in_line: point that is in the line
    # Output
    # Distance between point and line

    # return the distance
    return dist2line( point-point_in_line, line_vector )
end
#-----------------------------------------------------------------------------

# Computes angle using Al-Kashi
#-----------------------------------------------------------------------------
function angleAlKash( a::T1 , b::T2, c::T3 ) where { T1 <: Real, T2 <: Real, T3 <: Real }
    # Arguments
    # a,b,c: sides of the triangle to apply Al-Kashi
    # Output:
    # - Angle in degrees between a and b sides

    # Computes the cosine of the angle and then invert it into an angle in degrees
    return acosd( ( a*a + b*b - c*c)/( 2*a*b ) )
end
#-----------------------------------------------------------------------------

end
