using CharFuncPricing, Test
using TaylorSeries

@testset "SZCumulantsJoshi" begin
    κ = 4.0
    θ = 0.25
    σ = 1.0
    ρ = -0.5
    v0 = 0.1
    τ = 0.01
    params = SchobelZhuParams(v0, κ, θ, ρ, σ)
    t = Taylor1(Float64, 4)
    cft = evaluateLogCharFunc(Complex, params, t, τ)
    c1, c2, c4 = computeCumulants(params, τ)
    println(c1, " ", c2 / 2, " ", c4 / (2 * 3 * 4), " ", cft)
    @test isapprox(imag(cft.coeffs[2]), c1, atol = 1e-12)
    @test isapprox(-real(cft.coeffs[3]) * 2, c2, atol = 1e-12)
    @test isapprox(cft.coeffs[5] * 2 * 3 * 4, c4, atol = 1e-16)
end

@testset "SchobelZhuPaper" begin
    #reference values from "Stochastic Volatility with an Ornstein-Uhlenbeck Process: An Extension" (1998)
    κ = 4.0
    θ = 0.2
    σ = 0.1
    ρ = -0.5
    v0 = 0.2
    τ = 0.5
    spot = 100.0
    r = 0.0953
    forward = spot * exp(r * τ)
    df = exp(-r * τ)
    params = SchobelZhuParams(v0, κ, θ, ρ, σ)
    m = 128
    l = 6
    pricer =
        makeCosCharFuncPricer(Complex, Float64, Float64(MathConstants.pi), params, τ, m, l)
    strikes = [90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0]
    refPrices = [15.29, 11.50, 8.24, 5.60, 3.58, 2.16, 1.22]
    for (strike, refPrice) in zip(strikes, refPrices)
        price = priceEuropean(pricer, true, strike, forward, df)
        println(strike, " ", price)
        @test isapprox(refPrice, price, atol = 0.01)
    end
    θ = 0.0001
    v0 = 0.15
    params = SchobelZhuParams(v0, κ, θ, ρ, σ)
    pricer =
        makeCosCharFuncPricer(Complex, Float64, Float64(MathConstants.pi), params, τ, m, l)
    refPrices = [14.22, 9.60, 5.37, 2.15, 0.50, 0.06, 0.004]
    for (strike, refPrice) in zip(strikes, refPrices)
        price = priceEuropean(pricer, true, strike, forward, df)
        println(strike, " ", price)
        @test isapprox(refPrice, price, atol = 0.01)
    end

end
