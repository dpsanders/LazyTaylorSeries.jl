

const CoeffDict{N,T} = Dict{ SVector{N,Int}, T }

"""Lazy Taylor series with N variables and coefficients of type T
`degree` is the max degree of monomials allowed.
If `degree == -1` then the Tayor series is potentially infinite."""
immutable Taylor{N,T,F}
    f::F
    coeffs::CoeffDict{N,T}
    degree::Int   # maximum order of monomial
end


function Taylor(N, T, f::F, degree=-1) where {F}
    t = Taylor{N, T, F}(f, CoeffDict{N,T}(), degree)
    dummy = t[SVector(ntuple(_->0, Val{N}))]  # compile getindex by calculating first coefficient
    return t
end

degree(x::SVector) = sum(x)
degree(f::Taylor) = f.degree

function getindex(t::Taylor{N,T,F}, index::SVector) where {N,T,F}

    if degree(index) > t.degree
        return zero(T)
    end

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

function degree_sum(f, g)
    if f.degree == -1 || g.degree == -1
        return -1  # infinite
    end

    return  f.degree + g.degree
end

function degree_max(f, g)
    if f.degree == -1 || g.degree == -1
        return -1  # infinite
    end

    return max(f.degree, g.degree)
end



# these are memoized, but should look at performance without memoizing perhaps
+(f::Taylor{N,T}, g::Taylor{N,T}) where {N,T} = Taylor(N, T,
                                                (t, i) -> f[i] + g[i],
                                                degree_max(f, g) )

-(f::Taylor{N,T}, g::Taylor{N,T}) where {N,T} = Taylor(N, T,
                                                        (t, i) -> f[i] - g[i],
                                                        degree_max(f, g) )

#
# # formulas from Warwick Tucker, *Validated Numerics*
#
function *(f::Taylor{N,T}, g::Taylor{N,T}) where {N,T}

    return Taylor(N, T,
            (t, index) -> begin

                coeff = zero(T)

                tuples = generate_tuples(index)

                for i in tuples
                    coeff += f[i] * g[index - i]
                end

                coeff
            end,
            degree_sum(f, g)
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

"Generate tuples with a given maximum degree"
function generate_tuples(N, deg)

    t = zero(SVector{N,Int})

    tuples = [t]

    which = N

    @inbounds while degree(t) <= deg && which > 0

        if degree(t) == deg
            t = setindex(t, 0, which)
            which -= 1

            if which == 0
                break
            end

        else
            which = N
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

const variable_names = ["x", "y", "z", "w"]



function evaluate!(f::Taylor{N,T}) where {N,T}
    tuples = generate_tuples(N, degree(f))

    for t in tuples
        f[t]
    end
end


function Base.show(io::IO, f::Taylor{N,T}) where {N, T}

    if degree(f) == 0
        print(io, f[zero(SVector{N,Int})])
        return
    end

    # todo: order by degree
    for t in sort(collect(keys(f.coeffs)), lt=lexless)
        value = f[t]

        if value == 0
            continue
        end

        if (value == 1 && iszero(t))
            print(io, value)
        end

        if value != 1
            print(io, value)
        end



        for i in 1:N

            if t[i] == 1
                print(io, variable_names[i])
            elseif t[i] > 1
                print(io, variable_names[i], "^", t[i])
            end
        end

        print(io, " + ")
    end
end




export x, y, z, o

function variable(N, T, i)
    vec = SVector(ntuple(j->(j==i ? 1 : 0), Val{N}))
    return Taylor(N, T, (t,index)->(index==vec ? 1 : 0), 1)
end

function constant(N, T)
    return Taylor(N, T, (t,i)->(i==zero(SVector{N,T}) ? 1 : 0), 0)
end

function variables(N, T)
    vars = [variable(N, T, i) for i in 1:N]

    for var in vars
        evaluate!(var)
    end

    return (constant(N,T), vars ...)
end


export o, x, y, z
o, x, y, z = variables(3, Float64)
