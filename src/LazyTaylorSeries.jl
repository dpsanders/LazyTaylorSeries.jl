module LazyTaylorSeries

using StaticArrays


export LazyTaylor1, Taylor, tt
export variable, variables, constant, degree, evaluate!



import Base:
    +, -, *,
    exp,
    getindex


include("LazyTaylor1.jl")
include("taylorN.jl")


end
