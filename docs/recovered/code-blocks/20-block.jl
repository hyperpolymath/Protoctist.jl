using StaticArrays
@inline muladd3(a::SVector{3,Float32}, b::SVector{3,Float32}) = a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
