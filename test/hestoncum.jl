using CharFuncPricing, TaylorSeries

@testset "CumulantsJoshi" begin
    κ = 4.0
    θ = 0.25
    σ = 1.0
    ρ = -0.5
    v0 = 0.01
    τ = 0.01
    params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
    cf = DefaultCharFunc{HestonParams{Float64},Taylor1{Complex}}(params)
    t = Taylor1(Float64, 4)
    cft = CharFuncPricing.evaluateLogCharFuncCui(cf, t, τ)
    c1, c2, c4 = computeCumulants(params, τ)
    println(c1, " ", c2 / 2, " ", c4 / (2 * 3 * 4), " ", cft)
    @test isapprox(imag(cft.coeffs[2]), c1, atol = 1e-12)
    @test isapprox(-real(cft.coeffs[3]) * 2, c2, atol = 1e-12)
    @test isapprox(cft.coeffs[5] * 2 * 3 * 4, c4, atol = 1e-16)
end
@testset "CumulantsAlan" begin
    κ = 4.0
    θ = 0.25
    σ = 1.0
    ρ = -0.5
    v0 = 0.04
    τ = 1.0
    params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
    cf = DefaultCharFunc{HestonParams{Float64},Taylor1{Complex}}(params)
    t = Taylor1(Float64, 4)
    cft = CharFuncPricing.evaluateLogCharFuncCui(cf, t, τ)
    c1, c2, c4 = computeCumulants(params, τ)
    println(c1, " ", c2 / 2, " ", c4 / (2 * 3 * 4), " ", cft)
    @test isapprox(imag(cft.coeffs[2]), c1, atol = 1e-12)
    @test isapprox(-real(cft.coeffs[3]) * 2, c2, atol = 1e-12)
    @test isapprox(cft.coeffs[5] * 2 * 3 * 4, c4, atol = 1e-15)
end
@testset "CumulantsNearZeroKappa" begin
    #Low accuracy with small kappa due to powers of kappa in the denominator.
    κ = 0.02
    θ = 0.25
    σ = 1.0
    ρ = -0.5
    v0 = 0.04
    τ = 1.0
    params = HestonParams{Float64}(v0, κ, θ, ρ, σ)
    cf = DefaultCharFunc{HestonParams{Float64},Taylor1{Complex}}(params)
    t = Taylor1(Float64, 4)
    cft = CharFuncPricing.evaluateLogCharFuncCui(cf, t, τ)
    c1, c2, c4 = computeCumulants(params, τ)
    println(c1, " ", c2 / 2, " ", c4 / (2 * 3 * 4), " ", cft)
    @test isapprox(imag(cft.coeffs[2]), c1, atol = 1e-12)
    @test isapprox(-real(cft.coeffs[3]) * 2, c2, atol = 1e-3)
    @test isapprox(cft.coeffs[5] * 2 * 3 * 4, c4, atol = 1e-3)
end

@testset "HestonCumulantsSmallKappa" begin
    #the closed form cancels catastrophically for small κτ, references from BigFloat Taylor expansion of the char. func.
    refC4 = Dict(1e-2 => 0.009974159883692182, 1e-3 => 0.01003791387894611, 1e-6 => 0.010045022404969867)
    for (κ, c4) in refC4
        cum = CharFuncPricing.computeCumulants(CharFuncPricing.HestonParams(0.025, κ, 0.025, -0.65, 0.4), 1.0)
        @test isapprox(cum[3], c4, rtol=1e-8)
    end
    p = CharFuncPricing.HestonParams(0.025, 1.0, 0.025, -0.65, 0.4)
    @test all(isapprox.(CharFuncPricing.computeCumulantsODE(p, 1.0), CharFuncPricing.computeCumulants(p, 1.0), rtol=1e-9))
    for κ in [0.0, 1e-10]
        cf = DefaultCharFunc(CharFuncPricing.HestonParams(0.025, κ, 0.025, -0.65, 0.4))
        pricer = makeCosCharFuncPricer(cf, 1.0, 256, 12)
        refPricer = ALCharFuncPricer(cf, n=64)
        for strike in [0.8, 1.0, 1.2]
            @test isapprox(priceEuropean(pricer, true, strike, 1.0, 1.0, 1.0), priceEuropean(refPricer, true, strike, 1.0, 1.0, 1.0), atol=1e-5)
        end
    end
end
