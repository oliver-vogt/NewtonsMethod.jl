using NewtonsMethod
using Test

# Define functions to be tested
m(x) = (x - 4)^3
m′(x) = 3*(x - 4)^2
h(x) = log(x)
h′(x) = 1/x
b(x) = 2*x + x^4
b′(x) = 2 + 4*x^3
g(x) = 2 + x^2
g′(x) = 2*x

@testset "NewtonsMethod.jl" begin
    # Run tests on newton root with analytical derivatives
    @test newtonroot(m, m′, x_0=0.5).value ≈ 4.0 rtol=1E-7
    @test newtonroot(h, h′, x_0=0.5).value ≈ 1.0 rtol=1E-7
    @test newtonroot(b, b′, x_0=0.5).value ≈ 0.0 rtol=1E-7
    @test newtonroot(b, b′, x_0=-2.5).value ≈ -2^(1/3) rtol=1E-7

    # Run tests on newton root with nummerical derivatives
    @test newtonroot(m, x_0=0.5).value ≈ 4.0 rtol=1E-7
    @test newtonroot(h, x_0=0.5).value ≈ 1.0 rtol=1E-7
    @test newtonroot(b, x_0=0.5).value ≈ 0.0 rtol=1E-7
    @test newtonroot(b, x_0=-2.5).value ≈ -2^(1/3) rtol=1E-7

    # Run tests on newton root with analytical derivatives, BigFloat
    @test newtonroot(m, m′, x_0=BigFloat(0.5)).value ≈ 4.0 rtol=1E-7
    @test newtonroot(h, h′, x_0=BigFloat(0.5)).value ≈ 1.0 rtol=1E-7
    @test newtonroot(b, b′, x_0=BigFloat(-2.5)).value ≈ -2^(1/3) rtol=1E-7

    # Run tests on newton root with nummerical derivatives, BigFloat
    @test newtonroot(m, x_0=BigFloat(0.5)).value ≈ 4.0 rtol=1E-7
    @test newtonroot(h, x_0=BigFloat(0.5)).value ≈ 1.0 rtol=1E-7
    @test newtonroot(b, x_0=BigFloat(-2.5)).value ≈ -2^(1/3) rtol=1E-7

    # Run tests with non-convergence function
    @test newtonroot(g,g′, x_0=0.5) == nothing
    @test newtonroot(g, x_0=0.5) == nothing

    # Run tests with changed maxiter and tol for some of the functions
    @test newtonroot(h, x_0=0.5, maxiter=5) == nothing
    @test newtonroot(b, b′, x_0=BigFloat(-2.5), maxiter=10_000).value ≈ -2^(1/3) 
    @test newtonroot(m, x_0=0.5, tol=10E-9).value ≈ 4.0 rtol=1E-7
    @test isapprox(newtonroot(m, m′, x_0=BigFloat(0.5), tol=10E-2).value,4.0)

end