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
