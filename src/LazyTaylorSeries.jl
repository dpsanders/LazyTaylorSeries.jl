module LazyTaylorSeries

using StaticArrays


export Taylor1, Taylor, tt, constant



import Base:
    +, -, *,
    exp,
    getindex


include("taylor1.jl")
include("taylorN.jl")


end
