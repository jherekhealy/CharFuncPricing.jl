using CharFuncPricing, Test, Printf
using StatsBase

@testset "Gauthier" begin
    τs = [1.0, 10.0]
    S0 = 61.90
    strikes = [0.7, 1.0, 1.3]
    r = 0.03
    params1 = HestonParams(0.6^2, 0.9, 0.1, -0.5, 0.1)
    params2 = HestonParams(0.7^2, 1.2, 0.15, -0.5, 0.2)
    cf = DefaultCharFunc(DoubleHestonParams(params1, params2))
    refPrices = [27.6092 19.4569 13.9299;
        45.2866 41.4006 38.2780
    ]
    for (i, τ) = enumerate(τs)
        df = exp(-r * τ)
        for (j, K) = enumerate(strikes)
            forward = S0 / df

            pricer = makeCosCharFuncPricer(cf, τ, relStrike=K / forward, tol=1e-6)
            price = priceEuropean(pricer, true, K * S0, forward, τ, df)
            pricer2 = JoshiYangCharFuncPricer(cf, τ, n=16)
            price2 = priceEuropean(pricer2, true, K * S0, forward, τ, df)
            pricer3 = FlinnCharFuncPricer(cf,τ, tTol = 1e-10, qTol=1e-10)
            price3 = priceEuropean(pricer3, true, K * S0, forward, τ, df)
            println(τ, " ", K, " ", price, " ", price2, " ", price3)
            @test isapprox(refPrices[i, j], price, atol=1e-2)
            @test isapprox(refPrices[i, j], price2, atol=1e-2)
            @test isapprox(refPrices[i, j], price3, atol=1e-2)
        end
    end
end