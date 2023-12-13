using Test, CharFuncPricing

@testset "CGMY-Fang-European" begin    
    params = CGMYParams(1.0,5.0,5.0,0.5)
    cf = DefaultCharFunc(params)
    pricerCos = CharFuncPricing.makeCosCharFuncPricer(cf,1.0,128,12)
    τ = 1.0
    df = exp(-0.1*τ)
    p =priceEuropean(pricerCos,true,100.0,100.0/df,τ,df)
    @test isapprox(19.812948843, p, atol=1e-7)
    params = CGMYParams(1.0,5.0,5.0,1.5)
    cf = DefaultCharFunc(params)
    pricerCos = CharFuncPricing.makeCosCharFuncPricer(cf,1.0,128,12)
    p =priceEuropean(pricerCos,true,100.0,100.0/df,τ,df)
    @test isapprox(49.790905469, p, atol=1e-7)
end
