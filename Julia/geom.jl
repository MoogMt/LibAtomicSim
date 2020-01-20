module geom

using LinearAlgebra

function distance( vec1::Vector{T1}, vec2::Vector{T2} ) where { T1 <: Real , T2 <: Real }
    return LinearAlgebra.norm(vec1-vec2)
end

function dist2Line( vector_2_line::Vector{T1}, line_vector::Vector{T2} ) where { T1 <: Real, T2 <: Real}
    return abs(dot(vector_2_line,line_vector))/LinearAlgebra.norm(line_vector)
end

function pointFromLine( point::Vector{T1}, line_vector::Vector{T2} , point_in_line::Vector{T3} ) where { T1 <: Real, T2 <: Real, T3 <: Real }
    return dist2line( point-point_in_line, line_vector )
end

function angleAlKash( a::T1 , b::T2, c::T3 ) where { T1 <: Real, T2 <: Real, T3 <: Real }
    return acosd((a*a+b*b-c*c)/(2*a*b))
end

end
