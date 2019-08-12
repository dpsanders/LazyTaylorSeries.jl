struct Taylor1{T,F,memoize}
    f::F
    coeffs::Vector{T}
end


Base.literal_pow(::typeof(^), t::Taylor1, n::Integer) = t^n

import Base: ^
^(t::Taylor1, n::Integer) = power_by_squaring(t, n)
# uses custom version of power_by_squaring in powers.jl file

# the function f must take *two* variables if it is memoized;
# the first is used as an explicit reference to the current object when necessary
# In principle this is independent of whether it is memoized


Taylor1(f::F, memoize) where {F} = Taylor1(Float64, f, memoize)

function Taylor1(T, f::F, memoize) where {F}
    t = Taylor1{T,F,Val{memoize}}(f, T[])
    dummy = t[0]  # compile getindex by calculating first coefficient
    return t
end

Taylor1(f::Function) = Taylor1(f, true)  # memoize by default

# this version of getindex is for non-memoized
getindex(t::Taylor1{T,F,Val{false}}, i::Int) where {T, F} = (t.f)(i)

# Memoized; use NaN to indicate value not yet calculated
function getindex(t::Taylor1{T,F,Val{true}}, i::Int) where {T, F}
    j = i + 1
    coeffs = t.coeffs

    @inbounds if j <= length(coeffs)
        if isnan(coeffs[j])
            coeffs[j] = (t.f)(t, i)  # pass in the object as the first argument to the function for those functions that are recursive
        end

        return coeffs[j]

    else  # too short
        current_length = length(coeffs)
        resize!(coeffs, j)
        @inbounds coeffs[current_length+1:end] .= NaN
        @inbounds coeffs[end] = (t.f)(t, i)

        @inbounds return coeffs[end]
    end
end

# tt is the independent variable; non-memoized (nothing stored in memory)
# tt = Taylor1( i::Int -> (i == 1) * 1.0, false )

# tt = Taylor1(Interval{Float64}, (t, i) -> (i == 1) ? (1..1) : (0..0), true)
# use as constant(3); also non-memoized
constant(c::Real) = Taylor1( i::Int -> (i == 0) ? c : 0.0, false )
constant(T, c::Real) = Taylor1(T, i::Int -> (i == 0) ? c : 0.0, false )


# these are memoized, but should look at performance without memoizing perhaps

# use promotion!

+(f::Taylor1{T}, g::Taylor1{T}) where {T} = Taylor1(T, (t, i) -> f[i] + g[i], true )
-(f::Taylor1{T}, g::Taylor1{T}) where {T} = Taylor1(T, (t, i) -> f[i] - g[i], true )

-(f::Taylor1{T}) where {T} = Taylor1(T, (t, i) -> -f[i], true )

-(a::Real, f::Taylor1{T}) where {T} = Taylor1(T, (t, i) -> (i == 0) ? a-f[0] : -f[i], true )
+(a::Real, f::Taylor1{T}) where {T} = Taylor1(T, (t, i) -> (i == 0) ? a+f[0] : +f[i], true )

+(f::Taylor1, a::Real) = a + f
-(f::Taylor1, a::Real) = f + (-a)

# formulas from Warwick Tucker, *Validated Numerics*

*(f::Taylor1{T}, g::Taylor1{T}) where {T} = Taylor1(T, (t, k) -> sum(f[i] * g[k-i] for i in 0:k), true)

*(a::Real, f::Taylor1{T}) where {T} = Taylor1(T, (t, i) -> a*f[i], true)
*(f::Taylor1, a::Real) = a * f

# self is a reference to the object exp(g), that is used recursively
function exp(g::Taylor1{T}) where {T}
    function f(self, k)
        k == 0 && return exp(g[0])
        # dummy = g[k]  # preallocate g
        return sum(i * g[i] * self[k-i] for i in 1:k) / k
    end

    return Taylor1(T, f, true)

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


## tt = Taylor1(Interval{Float64}, i -> (i == 1) ? (1..1) : (0..0), false)
