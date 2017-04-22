module LazyTaylorSeries

export Taylor, tt, constant

import Base:
    +, -, *,
    exp,
    resize!,
    getindex

immutable Taylor{T,F,memoize}
    f::F
    coeffs::Vector{T}
end

# the function f must take *two* variables if it is memoized;
# the first is used as an explicit reference to the current object when necessary
# In principle this is independent of whether it is memoized

function Taylor{F}(f::F, memoize)
    t = Taylor{Float64,F,Val{memoize}}(f, Float64[])
    dummy = t[0]  # compile getindex by calculating first coefficient
    return t
end

Taylor(f::Function) = Taylor(f, true)  # memoize by default

# this version of getindex is for non-memoized
getindex{T,F}(t::Taylor{T,F,Val{false}}, i::Int) = (t.f)(i)

"""Resize a Taylor object to fit coefficients up to order n"""
resize!{T,F}(t::Taylor{T,F,Val{false}}, n::Int) = return

function resize!{T,F}(t::Taylor{T,F,Val{true}}, n::Int)
    new_size = n+1
    current_length = length(t.coeffs)

    new_size <= current_length && return

    resize!(t.coeffs, new_size)
    @inbounds t.coeffs[current_length+1:end] .= NaN
end

# Memoized; use NaN to indicate value not yet calculated
function getindex{T,F}(t::Taylor{T,F,Val{true}}, i::Int)
    j = i + 1
    coeffs = t.coeffs

    resize!(t, i)

    if isnan(coeffs[j])
        coeffs[j] = (t.f)(t, i)  # pass in the object as the first argument to the function for those functions that are recursive

        #@show object_id(t), i
    end

    return coeffs[j]
end

# tt is the independent variable; non-memoized (nothing stored in memory)
tt = Taylor( i::Int -> (i == 1) * 1.0, false )

# use as constant(3); also non-memoized
constant(c::Float64) = Taylor( i::Int -> (i == 0) ? c : 0.0, false )
constant(c::Real) = constant(Float64(c))


# these are memoized, but should look at performance without memoizing perhaps
+(f::Taylor, g::Taylor) = Taylor( (t, i) -> f[i] + g[i], true )
-(f::Taylor, g::Taylor) = Taylor( (t, i) -> f[i] - g[i], true )


# formulas from Warwick Tucker, *Validated Numerics*

*(f::Taylor, g::Taylor) = Taylor( (t, k) -> sum(f[i] * g[k-i] for i in 0:k), true)

# self is a reference to the object exp(g), that is used recursively
function exp(g::Taylor)
    function f(self, k)
        k == 0 && return exp(g[0])
        # dummy = g[k]  # preallocate g

        resize!(self, k)
        resize!(g, k)

        return sum(i * g[i] * self[k-i] for i in 1:k) / k
    end

    return Taylor(f)

end

end
