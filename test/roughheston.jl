using CharFuncPricing, Test
import SpecialFunctions: gamma

@testset "RoughHeston-HestonLimit" begin
    #H=1/2 rough Heston (no mean reversion) is Heston with κ=0 and a flat forward variance curve.
    τ = 1.0
    strikes = [0.8, 1.0, 1.2]
    params = RoughHestonParams(0.5, -0.65, 0.4, s -> 0.025)
    hParams = CharFuncPricing.HestonParams(0.025, 1e-10, 0.025, -0.65, 0.4)
    refPricer = ALCharFuncPricer(DefaultCharFunc(hParams), n=64)
    refPrices = [priceEuropean(refPricer, true, k, 1.0, τ, 1.0) for k in strikes]
    cfa = AdamsCharFunc(params, 1000)
    pricer = makeCosCharFuncPricer(cfa, τ, 256, -6.0, 6.0)
    @test isapprox([priceEuropean(pricer, true, k, 1.0, τ, 1.0) for k in strikes], refPrices, atol=1e-4)
    #the Black control variate is meant for Fourier-integral (Lewis) pricers, not Cos
    pricer = AdaptiveFilonCharFuncPricer(RoughHestonCVCharFunc(cfa), τ, qTol=1e-8, myTrans=CharFuncPricing.IdentityTransformation(0.0, 160.0))
    @test isapprox([priceEuropean(pricer, true, k, 1.0, τ, 1.0) for k in strikes], refPrices, atol=1e-4)
end

@testset "RoughHeston-PadeAdams" begin
    τ = 1.0
    H = 0.05
    nu = 0.8 * sqrt(2H) * gamma(H + 0.5)
    params = RoughHestonParams(H, -0.65, nu, s -> 0.025)
    strikes = [0.8, 1.0, 1.2]
    pricerPade = makeCosCharFuncPricer(PadeCharFunc(params), τ, 256, -6.0, 6.0)
    pricerAdams = makeCosCharFuncPricer(AdamsCharFunc(params, 1000), τ, 256, -6.0, 6.0)
    pricesPade = [priceEuropean(pricerPade, true, k, 1.0, τ, 1.0) for k in strikes]
    pricesAdams = [priceEuropean(pricerAdams, true, k, 1.0, τ, 1.0) for k in strikes]
    @test isapprox(pricesPade, pricesAdams, atol=2e-4)
    pricerPade64 = makeCosCharFuncPricer(PadeCharFunc{Complex}(params, FejerQuadrature(64)), τ, 256, -6.0, 6.0)
    @test isapprox([priceEuropean(pricerPade64, true, k, 1.0, τ, 1.0) for k in strikes], pricesPade, atol=1e-5)
end

@testset "PricerFactory" begin
    cf = DefaultCharFunc(CharFuncPricing.HestonParams(0.04, 1.5, 0.04, -0.6, 0.5))
    τ = 1.0
    ref = priceEuropean(ALCharFuncPricer(cf, n=64), true, 1.0, 1.0, τ, 1.0)
    for method in [:AndersenLake, :JoshiYang, :Flinn, :FlinnTransformed, :Swift, :Cos, :Cos128, :CosJunike, :FCC, "Joshi-Yang"]
        pricer = makeCharFuncPricer(method, cf, τ)
        @test isapprox(priceEuropean(pricer, true, 1.0, 1.0, τ, 1.0), ref, atol=1e-4)
    end
end
