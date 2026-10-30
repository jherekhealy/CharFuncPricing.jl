using CharFuncPricing, Test

@testset "FCC" begin
    cases = [
        (HestonParams(0.04, 4.0, 0.25, -0.5, 1.0), 1.0, [80.0, 100.0, 120.0]),
        (HestonParams(0.07, 0.3, 0.1, -0.1, 0.6), 1 / 12, [80.0, 100.0, 120.0]),
        (HestonParams(0.04, 0.5, 0.04, -0.9, 1.0), 10.0, [50.0, 100.0, 200.0]),
        (HestonParams(0.04, 1.0, 0.04, -0.7, 0.5), 0.01, [95.0, 100.0, 105.0]),
        (HestonParams(0.0001, 2.0, 0.0001, -0.5, 0.01), 0.25, [99.0, 100.0, 101.0]),
    ]
    forward = 100.0
    df = 0.98
    for (params, τ, strikes) in cases
        cf = DefaultCharFunc(params)
        refPricer = ALCharFuncPricer(cf, n=400)
        for isCall in [false, true]
            refPrices = [priceEuropean(refPricer, isCall, k, forward, τ, df) for k in strikes]
            pricer = FCCCharFuncPricer(cf, τ, qTol=1e-12)
            @test isapprox(priceEuropean(pricer, fill(isCall, length(strikes)), strikes, forward, τ, df), refPrices, atol=1e-12)
            pricer = FCCCharFuncPricer(cf, τ, qTol=1e-8)
            @test isapprox([priceEuropean(pricer, isCall, k, forward, τ, df) for k in strikes], refPrices, atol=1e-8)
        end
    end
    cf = DefaultCharFunc(HestonParams(0.04, 4.0, 0.25, -0.5, 1.0))
    τ = 1.0
    refPricer = ALCharFuncPricer(cf, n=400)
    #Black-Scholes control variate
    pricer = FCCCharFuncPricer(HestonCVCharFunc(cf, τ), τ, qTol=1e-10)
    @test isapprox(priceEuropean(pricer, false, 90.0, forward, τ, df), priceEuropean(refPricer, false, 90.0, forward, τ, df), atol=1e-10)
    #explicit truncation
    pricer = FCCCharFuncPricer(cf, τ, b=200.0)
    @test pricer.b == 200.0
    @test isapprox(priceEuropean(pricer, true, 110.0, forward, τ, df), priceEuropean(refPricer, true, 110.0, forward, τ, df), atol=1e-9)
    @test_throws DomainError priceEuropean(pricer, true, 100.0, forward, 2τ, df)
    #tolerance below the round-off level of the integrand: the refinement must stop
    cfLow = DefaultCharFunc(HestonParams(0.0001, 2.0, 0.0001, -0.5, 0.01))
    for N in [12, 16, 24], hp in [true, false]
        pricer = FCCCharFuncPricer(cfLow, 0.25, qTol=1e-14, N=N, hp=hp)
        @test length(pricer.centers) < 100
    end
    #uniform initial partition
    pricer = FCCCharFuncPricer(cf, τ, qTol=1e-10, dyadic=false)
    @test isapprox(priceEuropean(pricer, false, 90.0, forward, τ, df), priceEuropean(refPricer, false, 90.0, forward, τ, df), atol=1e-10)
    #Variance Gamma
    cfvg = DefaultCharFunc(CharFuncPricing.VarianceGammaParams(0.12, -0.14, 0.2))
    pricerCos = makeCosCharFuncPricer(cfvg, τ, tol=1e-10)
    pricer = FCCCharFuncPricer(cfvg, τ, qTol=1e-10)
    for k in [80.0, 100.0, 120.0]
        @test isapprox(priceEuropean(pricer, false, k, forward, τ, df), priceEuropean(pricerCos, false, k, forward, τ, df), atol=1e-9)
    end
end

@testset "FilonFlinnLowStrikes" begin
    #strikes below the forward give negative Filon/Flinn angles, the small-angle expansion must use |θ|
    cf = DefaultCharFunc(HestonParams(0.04, 0.5, 0.04, -0.9, 1.0))
    τ = 10.0
    refPricer = ALCharFuncPricer(cf, n=400)
    strikes = [5.0, 20.0, 50.0]
    refPrices = [priceEuropean(refPricer, false, k, 100.0, τ, 1.0) for k in strikes]
    for pricer in [AdaptiveFlinnCharFuncPricer(cf, τ, qTol=1e-10), AdaptiveFilonCharFuncPricer(cf, τ, qTol=1e-10)]
        @test isapprox([priceEuropean(pricer, false, k, 100.0, τ, 1.0) for k in strikes], refPrices, atol=1e-5)
    end
end

@testset "FCCNonFinite" begin
    #the Adams scheme with few steps is not finite for large u: FCC must stop the truncation there, and not hang
    params = RoughHestonParams(0.1, -0.65, 0.3, s -> 0.025)
    cfa = AdamsCharFunc(params, 100)
    pricer = @test_logs (:warn,) match_mode = :any FCCCharFuncPricer(cfa, 1.0, qTol=1e-8)
    @test pricer.b <= 256
    pricePade = priceEuropean(FCCCharFuncPricer(PadeCharFunc(params), 1.0, qTol=1e-10), true, 1.0, 1.0, 1.0, 1.0)
    @test isapprox(priceEuropean(pricer, true, 1.0, 1.0, 1.0, 1.0), pricePade, atol=1e-3)
    #explicit truncation in the non-finite region
    @test_throws DomainError FCCCharFuncPricer(cfa, 1.0, b=4096.0)
end
