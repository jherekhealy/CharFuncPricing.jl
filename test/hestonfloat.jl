using CharFuncPricing, Test
using TaylorSeries

#Testing https://financepress.com/2019/02/15/heston-model-reference-prices/
@testset "AlanFloat64Set" begin
    r = 0.01
    q = 0.02
    κ = 4.0
    θ = 0.25
    σ = 1.0
    ρ = -0.5
    v0 = 0.04
    τ = 1.0
    spot = 100.0
    strikes = [80.0, 90.0, 100.0, 110.0, 120.0]
    alanPuts = [
        7.958878113256768285213263077598987193482161301733,
        12.017966707346304987709573290236471654992071308187,
        17.055270961270109413522653999411000974895436309183,
        23.017825898442800538908781834822560777763225722188,
        29.811026202682471843340682293165857439167301370697,
    ]
    alanCalls = [
        26.774758743998854221382195325726949201687074848341,
        20.933349000596710388139445766564068085476194042256,
        16.070154917028834278213466703938231827658768230714,
        12.132211516709844867860534767549426052805766831181,
        9.024913483457835636553375454092357136489051667150,
    ]
    spot *= exp((r - q) * τ)
    df = exp(-r * τ)
    params = HestonParams(v0, κ, θ, ρ, σ)
    l = 32
    m = 1024
    pricer = makeCosCharFuncPricer(Complex, Float64, Float64(MathConstants.pi), params, τ, m, l)
    for (strike, refCall, refPut) in zip(strikes, alanCalls, alanPuts)
        price = priceEuropean(pricer, false, strike, spot, df)
        println(Float64(strike), " P ", price, " ", price - refPut)
        @test isapprox(Float64(price - refPut), 0, atol = 1e-13)
        price = priceEuropean(pricer, true, strike, spot, df)
        println(Float64(strike), " C ", price, " ", price - refCall)
        @test isapprox(Float64(price - refCall), 0, atol = 1e-13)
    end
end
@testset "AlanJoshiFloat64Set" begin
    r = 0.01
    q = 0.02
    κ = 4.0
    θ = 0.25
    σ = 1.0
    ρ = -0.5
    v0 = 0.01
    τ = 0.01
    spot = 100.0
    strikes = [90.0, 95.0, 100.0, 105.0, 110.0]
    alanPuts = [
        4.5183603586861772614990106188215872180542e-8,
        0.000461954855653851579672612557018857858641926937,
        0.477781171629504680023239655436072890669645669297,
        5.009501052563650299130635110520904481889436667608,
        10.008998550115123724684210555728039829315964456261,
    ]
    alanCalls = [
        9.989001595065276544935948045293485530832966049263,
        4.989963479738160122154264702582719627807098780529,
        0.467782671512844263098248405184095087949465507760,
        2.527447823194706060519991248106500619490942e-6,
        1.29932760052624920704881258510264466e-13,
    ]
    spot *= exp((r - q) * τ)
    df = exp(-r * τ)
    params = HestonParams(v0, κ, θ, ρ, σ)
    l = 32
    m = 1024
    pricer = makeCosCharFuncPricer(Complex, Float64, Float64(MathConstants.pi), params, τ, m, l)
    for (strike, refCall, refPut) in zip(strikes, alanCalls, alanPuts)
        price = priceEuropean(pricer, false, strike, spot, df)
        println(Float64(strike), " P ", price, " ", price - refPut)
        @test isapprox(Float64(price - refPut), 0, atol = 2e-14)
        price = priceEuropean(pricer, true, strike, spot, df)
        println(Float64(strike), " C ", price, " ", price - refCall)
        @test isapprox(Float64(price - refCall), 0, atol = 2e-14)
    end
end
@testset "CumulantsJoshi" begin
    κ = 4.0
    θ = 0.25
    σ = 1.0
    ρ = -0.5
    v0 = 0.01
    τ = 0.01
    params = HestonParams(v0, κ, θ, ρ, σ)
    t = Taylor1(Float64, 4)
    cft = evaluateLogCharFunc(Complex, params, t, τ)
    c1, c2, c4 = computeCumulants(params, τ)
    println(c1, " ", c2 / 2, " ", c4 / (2 * 3 * 4), " ", cft)
    @test isapprox(imag(cft.coeffs[2]), c1, atol = 1e-12)
    @test isapprox(-real(cft.coeffs[3]) * 2, c2, atol = 1e-12)
    @test isapprox(cft.coeffs[5] * 2 * 3 * 4, c4, atol = 1e-16)
end
@testset "CumulantsAlan" begin
    κ = 4.0
    θ = 0.25
    σ = 1.0
    ρ = -0.5
    v0 = 0.04
    τ = 1.0
    params = HestonParams(v0, κ, θ, ρ, σ)
    t = Taylor1(Float64, 4)
    cft = evaluateLogCharFunc(Complex, params, t, τ)
    c1, c2, c4 = computeCumulants(params, τ)
    println(c1, " ", c2 / 2, " ", c4 / (2 * 3 * 4), " ", cft)
    @test isapprox(imag(cft.coeffs[2]), c1, atol = 1e-12)
    @test isapprox(-real(cft.coeffs[3]) * 2, c2, atol = 1e-12)
    @test isapprox(cft.coeffs[5] * 2 * 3 * 4, c4, atol = 1e-16)
end
@testset "CumulantsNearZeroKappa" begin
    #Low accuracy with small kappa due to powers of kappa in the denominator.
    κ = 0.02
    θ = 0.25
    σ = 1.0
    ρ = -0.5
    v0 = 0.04
    τ = 1.0
    params = HestonParams(v0, κ, θ, ρ, σ)
    t = Taylor1(Float64, 4)
    cft = evaluateLogCharFunc(Complex, params, t, τ)
    c1, c2, c4 = computeCumulants(params, τ)
    println(c1, " ", c2 / 2, " ", c4 / (2 * 3 * 4), " ", cft)
    @test isapprox(imag(cft.coeffs[2]), c1, atol = 1e-12)
    @test isapprox(-real(cft.coeffs[3]) * 2, c2, atol = 1e-3)
    @test isapprox(cft.coeffs[5] * 2 * 3 * 4, c4, atol = 1e-3)
end
