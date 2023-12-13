using CharFuncPricing, Test

@testset "HestonDigital" begin
    hestonParams = HestonParams(0.007569, 3.46, 0.008, -0.82, 0.14)
    cf = DefaultCharFunc(hestonParams)
    pricerCos = CharFuncPricing.makeCosCharFuncPricer(cf,1.0,1024,16)
    pal = ALCharFuncPricer(cf)
    τ = 1.0
    strikes = [0.7,1.0,1.3]
    forward = 1.0
    for strike = strikes
        p = priceDigital(pricerCos,false,strike,forward,τ,1.0)
        pd = (priceEuropean(pal,false, strike+1e-4, forward, τ,1.0)-priceEuropean(pal,false, strike-1e-4, forward, τ,1.0))/2e-4
        println(strike," ", p, " ", pd-p)
        @test isapprox(p, pd, atol=3e-7)
        p = priceDigital(pricerCos,true,strike,forward,τ,1.0)
        pd = -(priceEuropean(pal,true, strike+1e-4, forward, τ,1.0)-priceEuropean(pal, true, strike-1e-4, forward, τ,1.0))/2e-4
        @test isapprox(p, pd, atol=3e-7)
    end

end
 
