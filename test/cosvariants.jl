using CharFuncPricing, Test

@testset "CosVariants" begin
    τ = 1.0
    forward = 100.0
    df = 0.99
    cf = DefaultCharFunc(HestonParams(0.04, 4.0, 0.25, -0.5, 1.0))
    refPricer = ALCharFuncPricer(cf, n=200)
    strikes = [80.0, 100.0, 120.0]
    refPrices(isCall) = [priceEuropean(refPricer, isCall, k, forward, τ, df) for k in strikes]
    prices(p, isCall) = [priceEuropean(p, isCall, k, forward, τ, df) for k in strikes]
    @testset "Lipton" begin
        p = makeCosLiptonCharFuncPricer(cf, τ, 256, 12)
        @test isapprox(prices(p, false), refPrices(false), atol=1e-10)
        @test isapprox(prices(p, true), refPrices(true), atol=1e-10)
        @test_throws DomainError priceEuropean(p, true, 100.0, forward, 2τ, df)
    end
    @testset "GilPelaez" begin
        p = CharFuncPricing.makeCosGPCharFuncPricer(cf, τ, 256, 12)
        @test isapprox(prices(p, false), refPrices(false), atol=1e-10)
        @test isapprox(prices(p, true), refPrices(true), atol=1e-10)
        #strikes outside of the truncation range
        @test priceEuropean(p, false, 1e-8, forward, τ, df) == 0
        @test isapprox(priceEuropean(p, false, 1e8, forward, τ, df), df * (1e8 - forward))
        @test_throws DomainError priceEuropean(p, true, 100.0, forward, 2τ, df)
    end
    @testset "Alpha" begin
        for α in [0.5, 1.0]
            p = makeCosAlphaCharFuncPricer(cf, τ, 256, 12, α)
            @test isapprox(prices(p, false), refPrices(false), atol=1e-6)
            @test isapprox(prices(p, true), refPrices(true), atol=1e-10)
        end
        p = makeCosAlphaCharFuncPricer(cf, τ, 256, 12, 0.5)
        @test_throws DomainError priceEuropean(p, true, 100.0, forward, 2τ, df)
    end
    @testset "FFT" begin
        p = CosFFTCharFuncPricer(cf, τ, 256, 12)
        isCalls = [false, true, false]
        ref = [refPrices(false)[1], refPrices(true)[2], refPrices(false)[3]]
        for isAlternative in [true, false]
            @test isapprox(priceEuropean(p, isCalls, strikes, forward, τ, df; isAlternative=isAlternative, useNFFT=false), ref, atol=1e-12)
            @test isapprox(priceEuropean(p, isCalls, strikes, forward, τ, df; isAlternative=isAlternative, useNFFT=true), ref, atol=1e-7)
        end
        @test_throws DomainError priceEuropean(p, isCalls, strikes, forward, 2τ, df)
    end
    @testset "CosTolerance" begin
        p = makeCosCharFuncPricer(cf, τ, tol=1e-8)
        @test isapprox(prices(p, false), refPrices(false), atol=1e-10)
        #number of terms from the Lambert W bound
        p = makeCosCharFuncPricer(cf, τ, tol=1e-8, s=0)
        @test isapprox(prices(p, false), refPrices(false), atol=1e-7)
    end
    @testset "CosEdgeCases" begin
        p = makeCosCharFuncPricer(cf, τ, 256, 12)
        @test priceEuropean(p, false, 1e-8, forward, τ, df) == 0
        @test isapprox(priceEuropean(p, false, 1e8, forward, τ, df), df * (1e8 - forward))
        @test priceDigital(p, false, 1e-8, forward, τ, df) == 0
        @test priceDigital(p, false, 1e8, forward, τ, df) == df
        @test_throws DomainError priceEuropean(p, true, 100.0, forward, 2τ, df)
        @test_throws DomainError priceDigital(p, true, 100.0, forward, 2τ, df)
        #generic cumulants by Taylor expansion of the characteristic function
        generic = invoke(CharFuncPricing.computeCumulants, Tuple{CharFuncPricing.CharFunc,Any}, cf, τ; n=2)
        @test all(isapprox.(generic, computeCumulants(cf, τ)[1:2], rtol=1e-12))
        @test_throws DomainError invoke(CharFuncPricing.computeCumulants, Tuple{CharFuncPricing.CharFunc,Any}, cf, τ; n=3)
        #characteristic functions without cumulants
        @test_throws ArgumentError makeCosCharFuncPricer(PadeCharFunc(RoughHestonParams(0.1, -0.65, 0.3, s -> 0.025)), τ, 64, 12)
    end
    @testset "ChebyshevCharFunc" begin
        ccf = CharFuncPricing.ChebyshevCharFunc(cf, 200.0, 128, τ)
        p = FlinnCharFuncPricer(ccf, τ)
        @test isapprox(prices(p, false), refPrices(false), atol=1e-7)
        acf = AdamsCharFunc(RoughHestonParams(0.1, -0.65, 0.3, s -> 0.025), 50)
        @test CharFuncPricing.getControlVariance(CharFuncPricing.ChebyshevCharFunc(acf, 50.0, 16, τ), τ) == CharFuncPricing.getControlVariance(acf, τ)
    end
    @testset "JoshiYang" begin
        p = JoshiYangCharFuncPricer(cf, τ, n=64)
        @test isapprox(prices(p, false), refPrices(false), atol=1e-12)
        @test isapprox(prices(p, true), refPrices(true), atol=1e-12)
        @test_throws DomainError priceEuropean(p, true, 100.0, forward, 2τ, df)
    end
end
