module LazyTaylorSeries

export Taylor, tt, constant

import Base:
    +, -, *,
    exp,
    getindex

immutable Taylor{T,F,memoize}
    f::F
    coeffs::Vector{T}
end

# the function f must take *two* variables
# the first is used as an explicit reference to the current object when necessary

function Taylor{F}(f::F, memoize)
    t = Taylor{Float64,F,Val{memoize}}(f, Float64[])
    dummy = t[0]  # compile getindex and
    return t
end

Taylor(f::Function) = Taylor(f, true)

getindex{T,F}(t::Taylor{T,F,Val{false}}, i::Int) = (t.f)(i)

# Memoized; use NaN to indicate value not yet calculated
function getindex{T,F}(t::Taylor{T,F,Val{true}}, i::Int)
    j = i + 1
    coeffs = t.coeffs

    @inbounds if j <= length(coeffs)
        if isnan(coeffs[j])
            coeffs[j] = (t.f)(t, i)
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

tt = Taylor( i::Int -> (i == 1) * 1.0, false )


constant(c::Float64) = Taylor( i::Int -> (i == 0) ? c : 0.0, false )
constant(c::Real) = constant(Float64(c))


+(f::Taylor, g::Taylor) = Taylor( (t, i) -> f[i] + g[i], true )
-(f::Taylor, g::Taylor) = Taylor( (t, i) -> f[i] - g[i], true )

# formulas from Warwick Tucker, *Validated Numerics*

@inbounds *(f::Taylor, g::Taylor) = Taylor( (t, k) -> sum(f[i] * g[k-i] for i in 0:k), true)



# self is a reference to the object exp(g)
function exp(g::Taylor)
    function f(self, k)
        k == 0 && return exp(g[0])
        # dummy = g[k]  # preallocate g
        return sum(i * g[i] * self[k-i] for i in 1:k) / k
    end

    return Taylor(f)

end

end
