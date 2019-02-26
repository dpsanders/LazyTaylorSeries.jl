


# struct Taylor1{T,F,memoize}
#     f::F
#     coeffs::Vector{T}
# end

struct Taylor1{T,memoize}
    f::Function
    coeffs::Dict{Int,T}
end

Base.literal_pow(::typeof(^), t::Taylor1, n::Integer) = Base.power_by_squaring(t, n)

import Base: ^
^(t::Taylor1, n::Integer) = Base.power_by_squaring(t, n)

# the function f must take *two* variables if it is memoized;
# the first is used as an explicit reference to the current object when necessary
# In principle this is independent of whether it is memoized

# function Taylor1(f::F, memoize) where F
#     t = Taylor1{Float64,F,Val{memoize}}(f, Float64[])
#     dummy = t[0]  # compile getindex by calculating first coefficient
#     return t
# end

function Base.show(io::IO, t::Taylor1{T,Val{true}}) where {T}
    for (k, v) in t.coeffs
        if k == 0
            print(io, "$v + ")
        elseif k == 1
            print(io, "$v t + ")
        else
            print(io, "$v t^$k + ")
        end
    end
end




function Taylor1(T, f, memoize)
    t = Taylor1{T,Val{memoize}}(f, Dict{Int,T}())
    dummy = t[0]  # compile getindex by calculating first coefficient
    return t
end

Taylor1{T}(f, memoize) where {T} = Taylor1(T, f, memoize)
Taylor1{T}(f) where {T} = Taylor1{T}(f, true)

# this version of getindex is for non-memoized
Base.getindex(t::Taylor1{T,Val{false}}, i::Int) where {T} = (t.f)(t, i)

# Memoized; use NaN to indicate value not yet calculated
function Base.getindex(t::Taylor1{T,Val{true}}, i::Int) where {T}
    coeffs = t.coeffs

    if !haskey(coeffs, i)
        coeffs[i] = (t.f)(t, i)
    end

    return coeffs[i]
end

# tt is the independent variable; non-memoized (nothing stored in memory)
# tt = Taylor1( (t, i) -> (i == 1) * 1.0, false )

# use as constant(3); also non-memoized
constant(c) = Taylor1( (t, i) -> (i == 0) ? c : zero(c), false )



# these are memoized, but should look at performance without memoizing perhaps

# use promotion!

+(f::Taylor1{T}, g::Taylor1{T}) where {T} = Taylor1{T}( (t, i) -> f[i] + g[i] )
-(f::Taylor1{T}, g::Taylor1{T}) where {T} = Taylor1{T}( (t, i) -> f[i] - g[i] )

-(f::Taylor1{T}) where {T} = Taylor1{T}( (t, i) -> -f[i] )

-(a::Real, f::Taylor1{T}) where {T} = Taylor1{T}( (t, i) -> (i == 0) ? a-f[0] : -f[i] )
+(a::Real, f::Taylor1{T}) where {T} = Taylor1{T}( (t, i) -> (i == 0) ? a+f[0] : +f[i] )


# formulas from Warwick Tucker, *Validated Numerics*

*(f::Taylor1{T}, g::Taylor1{T}) where {T} = Taylor1( (t, k) -> sum(f[i] * g[k-i] for i in 0:k) )

*(a::Real, f::Taylor1{T}) where {T} = Taylor1( (t, i) -> a*f[i] )

# self is a reference to the object exp(g), that is used recursively
function exp(g::Taylor1{T}) where {T}
    function f(self, k)
        k == 0 && return exp(g[0])
        # dummy = g[k]  # preallocate g
        return sum(i * g[i] * self[k-i] for i in 1:k) / k
    end

    return Taylor1{T}(f)

end

"""
Evaluate using Horner rule
"""
function (f::Taylor1)(x)
    total = zero(x)

    for (k, v) in f.coeffs
        total += v * x^k
    end

    return total
end
