"""Lazy Taylor series in 1 variable.
`T` is the type of the coefficients.
`coeffs` are the stored (memoized) coefficients.
`f` is the function that calculates coefficients.

Note: `f` must take *two* variables;
the first is used as an explicit reference to the current object when necessary, e.g. for `exp`.
"""

using StaticArrays

const ContainerType = Union{Vector{T}, Dict{Int,T}, MVector{N,T}} where {T,N}


# similar to eltype:
coeff_type(::Type{Vector{T}}) where {T} = T
coeff_type(::Type{Dict{Int,T}}) where {T} = T

struct Taylor1{T,F,C<:ContainerType}
    f::F
    coeffs::C  # C is the type of the container
end

init(::Type{Vector{T}}) where {T} = T[]
init(::Type{Dict{Int,T}}) where {T} = Dict{Int,T}()

Taylor1(f::F, coeffs::C) where {F,C} = Taylor1{coeff_type(C),F,C}(f, coeffs)

Base.haskey(v::Vector, i::Int) = i in keys(v)   # type piracy

# Taylor1(T, f::F, C) where {F} = Taylor1{T,F,C}(f, init(C))

# Taylor1(f::F, C) where {F} = Taylor1(coeff_type(C), f, init(C))
# function Taylor1(f::F) where F
#     t = Taylor1{Float64,F}(f, Float64[])
#     dummy = t[0]  # compile getindex by calculating first coefficient
#     return t
# end

# function Taylor1(T, f::F) where {F}
#     t = Taylor1(f, T[])
#     #dummy = t[0]  # compile getindex by calculating first coefficient
#     return t
# end

# Taylor1(T, f::Function) = Taylor1(T, f)


# Memoized; use NaN to indicate value not yet calculated
# function getindex(t::Taylor1{T,F,C}, i::Int) where {T, F, C}
#     j = i + 1
#     coeffs = t.coeffs
#
#     @inbounds if j <= length(coeffs)
#         if isnan(coeffs[j])
#             coeffs[j] = (t.f)(t, i)  # pass in the object as the first argument to the function for those functions that are recursive
#         end
#
#         return coeffs[j]
#
#     else  # too short
#         current_length = length(coeffs)
#         resize!(coeffs, j)
#         @inbounds coeffs[current_length+1:end] .= NaN
#         @inbounds coeffs[end] = (t.f)(t, i)
#
#         @inbounds return coeffs[end]
#     end
# end


function Base.setindex!(t::Taylor1{T,F,Vector{T}}, val::T, i::Int) where {T,F}
    j = i + 1
    coeffs = t.coeffs

    current_length = length(coeffs)
    resize!(coeffs, j)
    coeffs[current_length+1:end] .= NaN
    coeffs[end] = val

    return coeffs[end]
end

function Base.setindex!(t::Taylor1{T,F,Dict{Int,T}}, val::T, i::Int) where {T,F}
    t.coeffs[i] = val
    return t.coeffs[i]
end


function Base.getindex(t::Taylor1{T,F,C}, i::Int) where {T, F, C}
    j = i + 1  # offset so that numbering is 0-based
    coeffs = t.coeffs

    if !haskey(coeffs, j)
        t[i] = (t.f)(t, i)
    end

    return coeffs[j]
end


# tt is the independent variable; non-memoized (nothing stored in memory)
# tt = Taylor1( i::Int -> (i == 1) * 1.0, false )

# tt = Taylor1(Interval{Float64}, (t, i) -> (i == 1) ? (1..1) : (0..0), true)
# use as constant(3); also non-memoized
constant(c::Real) = Taylor1( i::Int -> (i == 0) ? c : 0.0, false )
constant(T, c::Real) = Taylor1(T, i::Int -> (i == 0) ? c : 0.0, false )


# these are memoized, but should look at performance without memoizing perhaps

# use promotion!

import Base: +, -, *

+(f::Taylor1{T,F1,C}, g::Taylor1{T,F2,C}) where {T,F1,F2,C} = Taylor1( (t, i) -> f[i] + g[i], init(C))
-(f::Taylor1{T}, g::Taylor1{T}) where {T} = Taylor1(T, (t, i) -> f[i] - g[i])

-(f::Taylor1{T}) where {T} = Taylor1(T, (t, i) -> -f[i])

-(a::Real, f::Taylor1{T}) where {T} = Taylor1(T, (t, i) -> (i == 0) ? a-f[0] : -f[i])
+(a::Real, f::Taylor1{T}) where {T} = Taylor1(T, (t, i) -> (i == 0) ? a+f[0] : +f[i])

+(f::Taylor1, a::Real) = a + f
-(f::Taylor1, a::Real) = f + (-a)

# formulas from Warwick Tucker, *Validated Numerics*

*(f::Taylor1{T,F1,C}, g::Taylor1{T,F2,C}) where {T,F1,F2,C} = Taylor1( (t, i) -> sum(f[k] * g[i-k] for k in 0:i), init(C))

*(a::Real, f::Taylor1{T}) where {T} = Taylor1(T, (t, i) -> a*f[i])
*(f::Taylor1, a::Real) = a * f

# self is a reference to the object exp(g), that is used recursively
function exp(g::Taylor1{T}) where {T}
    function f(self, k)
        k == 0 && return exp(g[0])
        # dummy = g[k]  # preallocate g
        return sum(i * g[i] * self[k-i] for i in 1:k) / k
    end

    return Taylor1(T, f)

end

"""
Evaluate using Horner rule
"""
function (f::Taylor1)(x)
    total = f.coeffs[end]

    for i in length(f.coeffs)-1 : -1 : 1
        total = x * total + f.coeffs[i]
    end

    return total
end


Base.literal_pow(::typeof(^), t::Taylor1, n::Integer) = t^n

import Base: ^
^(t::Taylor1, n::Integer) = power_by_squaring(t, n)  # uses modified version in powers.jl
