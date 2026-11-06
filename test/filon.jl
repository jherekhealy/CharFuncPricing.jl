using Test, CharFuncPricing

@testset "Filon-truncation" begin
    #Filon with the truncation of each model (a bound on the frequency) against the adaptive Filon
    cases = [(DefaultCharFunc(HestonParams(0.04, 1.5, 0.04, -0.7, 0.6)), 1.0),
        (DefaultCharFunc(VarianceGammaParams(-0.14, 0.12, 0.2)), 1.0),
        (DefaultCharFunc(CGMYParams(1.0, 5.0, 5.0, 0.5)), 1.0),
        (DefaultCharFunc(CGMYParams(1.0, 5.0, 5.0, 1.5)), 1.0),
        (DefaultCharFunc(CGMYParams(0.5, 2.0, 3.5, 1.1, 0.1)), 0.25),
        (DefaultCharFunc(KouParams(0.16, 1.0, 0.4, 10.0, 5.0)), 0.5),
        (DefaultCharFunc(HestonKouParams(HestonParams(0.04, 1.5, 0.04, -0.7, 0.6), 0.5, 0.3, 25.0, 10.0)), 1.0)]
    for (cf, τ) in cases
        pf = FilonCharFuncPricer(cf, τ)
        pa = AdaptiveFilonCharFuncPricer(cf, τ, qTol=1e-12)
        for K in (0.7, 0.9, 1.0, 1.1, 1.4)
            @test isapprox(priceEuropean(pf, true, K, 1.0, τ, 1.0), priceEuropean(pa, true, K, 1.0, τ, 1.0), atol=1e-7)
        end
    end
end
