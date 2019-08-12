module LazyTaylorSeries

using StaticArrays


export Taylor1, Taylor, tt
export variable, variables, constant, degree, evaluate!



import Base:
    +, -, *,
    exp,
    getindex

include("powers.jl")
include("taylor1.jl")
include("taylorN.jl")


end
