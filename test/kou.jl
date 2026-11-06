using Test, CharFuncPricing
import CharFuncPricing: evaluateLogCharFunc, computeCumulants

@testset "Kou-Kou2002" begin
    #example of Kou (2002): S=100, K=98, T=0.5, r=5%, σ=16%, λ=1, p=0.4, η₁=10, η₂=5, call price 9.14732
    cf = DefaultCharFunc(KouParams(0.16, 1.0, 0.4, 10.0, 5.0))
    τ = 0.5; F = 100 * exp(0.05τ); df = exp(-0.05τ)
    @test isapprox(priceEuropean(AdaptiveFilonCharFuncPricer(cf, τ, qTol=1e-12), true, 98.0, F, τ, df), 9.14732, atol=1e-5)
    @test isapprox(priceEuropean(makeCosCharFuncPricer(cf, τ, 512, 16), true, 98.0, F, τ, df), 9.14732, atol=1e-5)
    @test isapprox(priceEuropean(makeSwiftCharFuncPricer(cf, τ, 8, 12), true, 98.0, F, τ, df), 9.14732, atol=1e-5)
    #closed-form cumulants against the Taylor expansion of the log characteristic function
    l = evaluateLogCharFunc(cf, CharFuncPricing.TaylorSeries.Taylor1(Float64, 5), τ)
    @test isapprox(collect(computeCumulants(cf, τ)), [imag(l.coeffs[2]), -2real(l.coeffs[3]), 24real(l.coeffs[5])], rtol=1e-12)
    #no jumps: Black-Scholes
    cf0 = DefaultCharFunc(KouParams(0.16, 0.0, 0.4, 10.0, 5.0))
    @test isapprox(priceEuropean(AdaptiveFilonCharFuncPricer(cf0, τ, qTol=1e-12), true, 98.0, F, τ, df), 6.9682847, atol=1e-6)
    @test_throws DomainError KouParams(0.16, 1.0, 0.4, 1.0, 5.0)
end

@testset "Kou-HestonKou" begin
    h = HestonParams(0.04, 1.5, 0.04, -0.7, 0.6)
    ph = AdaptiveFilonCharFuncPricer(DefaultCharFunc(h), 1.0, qTol=1e-12)
    p0 = AdaptiveFilonCharFuncPricer(DefaultCharFunc(HestonKouParams(h, 0.0, 0.3, 25.0, 10.0)), 1.0, qTol=1e-12)
    @test isapprox(priceEuropean(p0, true, 0.9, 1.0, 1.0, 1.0), priceEuropean(ph, true, 0.9, 1.0, 1.0, 1.0), atol=1e-14)
    cf = DefaultCharFunc(HestonKouParams(h, 0.5, 0.3, 25.0, 10.0))
    pf = AdaptiveFilonCharFuncPricer(cf, 1.0, qTol=1e-13); pc = makeCosCharFuncPricer(cf, 1.0, 1024, 16)
    for K in (0.6, 1.0, 1.4)
        @test isapprox(priceEuropean(pf, true, K, 1.0, 1.0, 1.0), priceEuropean(pc, true, K, 1.0, 1.0, 1.0), atol=1e-10)
    end
    #the cumulants add: Heston and the jumps of Kou
    c = computeCumulants(cf, 1.0); ch = computeCumulants(h, 1.0)
    @test c[2] > ch[2] && c[3] > ch[3]
end
