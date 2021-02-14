using CharFuncPricing, Test
using TaylorSeries

@testset "SZLogContinuity" begin
 v0 = 0.08
 ρ=-0.8
 σ=0.25
 θ=0.1
 κ=3.0
 #params from Cui et al for Heston.
 params = SchobelZhuParams{Float64}(sqrt(v0), κ, sqrt(θ), ρ, σ)
 cf = DefaultCharFunc(params)
 #with the original formula from SZ, there is a discontinuity in the log.
 sz = CharFuncPricing.evaluateLogCharFunc( cf, 2.0+0.0im, 15.0)
 lk = CharFuncPricing.evaluateLogCharFuncLK( cf, 2.0+0.0im, 15.0)
 @test isapprox(lk,sz, atol=1e-13)
 for u = 0.1:0.1:4.0
      lk = CharFuncPricing.evaluateLogCharFuncLK(cf, u+0.0im, 15.0)
      sz=CharFuncPricing.evaluateLogCharFunc( cf, u+0.0im, 15.0)
      @test isapprox(lk,sz, atol=1e-13)
  end

  v0 = 0.010201
  ρ  = -0.7
  σ = 0.61
  κ = 6.21
  θ = 0.019
  α = 3.35861
  params = SchobelZhuParams{Float64}(sqrt(v0), κ, sqrt(θ), ρ, σ)
  cf = DefaultCharFunc(params)
   #with the original formula from SZ, there is a discontinuity in the log.
  sz = CharFuncPricing.evaluateLogCharFunc(cf, 10.0+α*1im, 10.0)
  lk = CharFuncPricing.evaluateLogCharFuncLK(cf, 10.0+α*1im, 10.0)
  @test isapprox(lk,sz, atol=1e-13)

end

@testset "SZCumulantsJoshi" begin
    κ = 4.0
    θ = 0.25
    σ = 1.0
    ρ = -0.5
    v0 = 0.1
    τ = 0.01
    params = SchobelZhuParams{Float64}(v0, κ, θ, ρ, σ)
    cf = DefaultCharFunc{SchobelZhuParams{Float64}, Taylor1{Complex}}(params)
    t = Taylor1(Float64, 4)
    cft = evaluateLogCharFunc(cf, t, τ)
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
    params = SchobelZhuParams{Float64}(v0, κ, θ, ρ, σ)
    cf = DefaultCharFunc(params)
    m = 128
    l = 6
    pricer =
        makeCosCharFuncPricer(cf, τ, m, l)
    strikes = [90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0]
    refPrices = [15.29, 11.50, 8.24, 5.60, 3.58, 2.16, 1.22]
    for (strike, refPrice) in zip(strikes, refPrices)
        price = priceEuropean(pricer, true, strike, forward, τ, df)
        println(strike, " ", price)
        @test isapprox(refPrice, price, atol = 0.01)
    end
    θ = 0.0001
    v0 = 0.15
    params = SchobelZhuParams{Float64}(v0, κ, θ, ρ, σ)
    cf = DefaultCharFunc(params)

    pricer =
        makeCosCharFuncPricer(cf, τ, m, l)
    refPrices = [14.22, 9.60, 5.37, 2.15, 0.50, 0.06, 0.004]
    for (strike, refPrice) in zip(strikes, refPrices)
        price = priceEuropean(pricer, true, strike, forward,τ, df)
        println(strike, " ", price)
        @test isapprox(refPrice, price, atol = 0.01)
    end

end
