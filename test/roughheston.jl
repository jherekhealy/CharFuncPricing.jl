using CharFuncPricing, Test, ForwardDiff
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

@testset "AdamsDuals" begin
    #ForwardDiff duals in u, as used by the Joshi-Yang control variate, and broadcasting over the char. func.
    params = RoughHestonParams(0.1, -0.65, 0.3, s -> 0.025)
    cf = AdamsCharFunc(params, 100)
    @test isfinite(ForwardDiff.derivative(u -> real(CharFuncPricing.evaluateCharFunc(cf, u - 1im, 1.0)), 0.3))
    cfp = PadeCharFunc(params)
    @test isapprox(priceEuropean(JoshiYangCharFuncPricer(cfp, 1.0, n=64), true, 1.0, 1.0, 1.0, 1.0),
        priceEuropean(FCCCharFuncPricer(cfp, 1.0, qTol=1e-10), true, 1.0, 1.0, 1.0, 1.0), atol=1e-8)
end

@testset "AdamsConvergence" begin
    #order 1 + α with the trapezoidal rule, Richardson extrapolation with this order
    params = RoughHestonParams(0.1, -0.65, 0.3, s -> 0.025)
    L(cf) = CharFuncPricing.evaluateLogCharFunc(cf, 1.0 - 0.5im, 1.0)
    r1, r2, r3 = L(AdamsCharFunc(params, 6400)), L(AdamsCharFunc(params, 3200)), L(AdamsCharFunc(params, 1600))
    order = log2(abs(r3 - r2) / abs(r2 - r1))
    @test isapprox(order, 1.6, atol=0.05)
    ref = r1 + (r1 - r2) / (2^order - 1)
    @test abs(L(AdamsCharFunc(params, 100)) - ref) < 1e-6
    @test abs(L(AdamsCharFunc(params, 100, richardson=true)) - ref) < abs(L(AdamsCharFunc(params, 200)) - ref) / 10
end

@testset "JoshiYangNonFinite" begin
    #the Adams scheme is not finite at the largest Gauss-Laguerre nodes, where the integrand is negligible
    params = RoughHestonParams(0.1, -0.65, 0.3, s -> 0.025)
    cfa = AdamsCharFunc(params, 200)
    pricer = JoshiYangCharFuncPricer(cfa, 1.0, n=64)
    ref = priceEuropean(FCCCharFuncPricer(cfa, 1.0, qTol=1e-10), true, 1.0, 1.0, 1.0, 1.0)
    @test isapprox(priceEuropean(pricer, true, 1.0, 1.0, 1.0, 1.0), ref, atol=1e-6)
    #a non finite value inside the integration domain is an error
    phi = ComplexF64[1.0, NaN, 1e-20]
    @test_throws DomainError CharFuncPricing.zeroNonFiniteTail!(phi, [1.0, 2.0, 3.0], [1.0, 1.0, 1.0])
end
