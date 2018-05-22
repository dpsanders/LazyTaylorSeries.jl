

const CoeffDict{N,T} = Dict{ SVector{N,Int}, T }

immutable Taylor{N,T,F}
    f::F
    coeffs::CoeffDict{N,T}
end


function Taylor(N, T, f::F) where {F}
    t = Taylor{N, Float64, F}(f, CoeffDict{N,T}())
    dummy = t[SVector(ntuple(_->0, Val{N}))]  # compile getindex by calculating first coefficient
    return t
end


function getindex(t::Taylor{N,T,F}, index::SVector) where {N,T,F}

    coeffs = t.coeffs

    if haskey(coeffs, index)
        return coeffs[index]
    end

    coeffs[index] = (t.f)(t, index)

    return coeffs[index]

end

function getindex(t::Taylor{N,T,F}, index...) where {N,T,F}
    getindex(t, SVector(index))
end


# these are memoized, but should look at performance without memoizing perhaps
+(f::Taylor{N,T}, g::Taylor{N,T}) where {N,T} = Taylor(N,T, (t, i) -> f[i] + g[i] )
-(f::Taylor{N,T}, g::Taylor{N,T}) where {N,T} = Taylor(N,T, (t, i) -> f[i] - g[i] )

#
# # formulas from Warwick Tucker, *Validated Numerics*
#
function *(f::Taylor{N,T}, g::Taylor{N,T}) where {N,T}

    return Taylor(N, T,
            (t, index) -> begin

                coeff = zero(T)

                tuples = generate_tuples(index)

                for t in tuples
                    coeff += f[t] * g[index - t]
                end

                coeff
            end
            )
end



    #     for monomial_f in keys(f.coeffs)
    #         for monomial_g in keys(g.coeffs)
    #
    #             monomial_product = monomial_f + monomial_g
    #
    #
    #


# enumerate tuples less than or equal to t, varying over dimension which

"Generate tuples that are lexicographically less than x"
function generate_tuples(x::SVector{N,Int}) where N

    t = zero(x)

    tuples = [t]

    @inbounds while t != x
        #@show t
        which = N

        while t[which] == x[which] && which > 1
            t = setindex(t, 0, which)
            which -= 1
        end

        t = setindex(t, t[which]+1, which)
        push!(tuples, t)
    end

    return tuples

end
#
# struct TupleGenerator{N,T}
#     x::SVector{N,T}
# end
#
# function Base.start(g::TupleGenerator{N,T}) where {N,T}
#     return zero(g.x)
# end
#
# function Base.next(g::TupleGenerator{N,T}, state) where {N,T}
#     current = state
#
#     which = N
#
#     while state[which] == g.x[which] && which > 1
#         state = setindex(state, 0, which)
#         which -= 1
#     end
#
#     state = setindex(state, state[which]+1, which)
#
#     return current, state
# end
#
# function Base.done(g::TupleGenerator{N,T}, state) where {N,T}
#     return state == g.x
# end

# enum_tuples(SVector(1,1), 1)



# self is a reference to the object exp(g), that is used recursively
# function exp(g::Taylor)
#     function f(self, k)
#         k == 0 && return exp(g[0])
#         # dummy = g[k]  # preallocate g
#         return sum(i * g[i] * self[k-i] for i in 1:k) / k
#     end
#
#     return Taylor(f)
#
# end
#
# end


 x = Taylor(2, Float64, (t,i)->(i==SVector(1, 0) ? 1 : 0))
 y = Taylor(2, Float64, (t,i)->(i==SVector(0, 1) ? 1 : 0))
