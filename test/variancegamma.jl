using Test, CharFuncPricing

@testset "Hirsa" begin
    S0=100.0
    σ = 0.12
    θ = -0.14
    ν = 0.2
    strikes = [90.0, 100.0, 120.0]
    cf = DefaultCharFunc(CharFuncPricing.VarianceGammaParams(θ, σ, ν))
    τ = 1.0
    r=0.1
    refPrices = [19.0994, 11.3700,1.92123]
    pricerCosFFT = CharFuncPricing.CosFFTCharFuncPricer(cf, τ, 1024, 8)
    prices = priceEuropean(pricerCosFFT, ones(Bool, length(strikes)), strikes, S0*exp(r*τ), τ, exp(-r*τ))
    @test isapprox(refPrices, prices, atol=2e-4)
    refPrice100 = 11.3700278104 #from crisostomo paper
    pricerCos = makeCosCharFuncPricer(cf, τ,256*4,8)
    price100 = priceEuropean(pricerCos, true, S0, S0 * exp(r*τ), τ, exp(-r*τ))
    @test isapprox(refPrice100, price100, atol=2e-9)

end
@testset "VG-AndersenLake-Case1And2" begin
    S0=100.0

    σ = 0.12136
    θ = -0.1436
    ν = 0.3
    r = 0.1
    strikes = [60.0, 101.0, 140.0]
    Ts = [0.1, 1.0]
    refPrices = [40.5972193355 1.3938439616 0.0000061410;
    45.7164396686 10.9815614276 0.101970645
        ]
    cf = DefaultCharFunc(CharFuncPricing.VarianceGammaParams(θ, σ, ν))
    prices = zeros(length(Ts), length(strikes))
    for (i, τ) in enumerate(Ts)
        pricerCosFFT = CharFuncPricing.CosFFTCharFuncPricer(cf, τ, 1024, 10)
        prices[i,:] = priceEuropean(pricerCosFFT, ones(Bool, length(strikes)), strikes, S0*exp(r*τ), τ, exp(-r*τ))
    end
    #what is the effective forward? AndersenLake modify the original Crisotomo test by using F=100 intead of S0=100. very similar to using r=0.0.

end

@testset "VG-AndersenLake-Case4And5" begin
    S0=100.0

    θ = 1.5
    σ = 1.0
    ν = 0.2
    r = 0.02
    strikes = [60.0, 90.0, 140.0]
    Ts = [0.1, 1.0]
    #refPrices = [40.5972193355 1.3938439616 0.0000061410;
    #45.7164396686 10.9815614276 0.101970645        ]
    cf = DefaultCharFunc(CharFuncPricing.VarianceGammaParams(θ, σ, ν))
    prices = zeros(length(Ts), length(strikes))
    for (i, τ) in enumerate(Ts)
        pricerCosFFT = CharFuncPricing.CosFFTCharFuncPricer(cf, τ, 1024, 10)
        prices[i,:] = priceEuropean(pricerCosFFT, ones(Bool, length(strikes)), strikes, S0*exp(r*τ), τ, exp(-r*τ))
    end
    

end
