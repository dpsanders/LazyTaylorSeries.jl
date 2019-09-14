module LazyTaylorSeries

using StaticArrays


export Taylor1, Taylor, tt
export variable, variables, constant, degree, evaluate!
export reset!


import Base:
    +, -, *,
    exp,
    getindex

include("powers.jl")
include("taylor1.jl")
include("taylorN.jl")

variable(T) = Taylor1(T, (t, i) -> (i == 1) ? one(T) : zero(T))


end
