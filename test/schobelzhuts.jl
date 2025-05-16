using CharFuncPricing, Test

@testset "SchobelZhuTSConstantPaper" begin
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
    params = SchobelZhuTSParams(v0, [κ], [θ], [ρ], [σ], [0.0])
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
    params = SchobelZhuTSParams(v0, [κ], [θ], [ρ], [σ], [0.0])
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
