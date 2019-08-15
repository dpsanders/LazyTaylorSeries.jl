using LazyTaylorSeries
using Test


@testset "Basic usage" begin
    t = Taylor1(Float64, (t, i) -> (i == 1))  # define variable
    t2 = variable(Float64)

    @test t[0] == 0 == t2[0]
    @test t[1] == 1 == t2[1]
    @test t[2] == 0 == t2[2]
    @test t[100] == 0 == t2[100]

    s = 1 + t
    @test (s[0], s[1], s[2], s[3]) == (1, 1, 0, 0)

    s = t + t
    @test (s[0], s[1], s[2], s[3]) == (0, 2, 0, 0)

    s = t * t
    @test (s[0], s[1], s[2], s[3]) == (0, 0, 1, 0)

    s = t^2
    @test (s[0], s[1], s[2], s[3]) == (0, 0, 1, 0)

    s = (1 + t)^2
    @test (s[0], s[1], s[2], s[3]) == (1, 2, 1, 0)

end
