using CharFuncPricing, Test
using TaylorSeries
using StatsBase
import CharFuncPricing: lambertW

 @testset "LF2y" begin
    r = 0.0
    q = 0.0
    κ = 1.0
    θ = 0.1
    σ = 1.0
    ρ = -0.9
    v0 = 0.1
    τ = 2.0
    spot = 1.0

    # κ = 1.417
    # θ = 0.082
    # σ = 0.591
    # ρ = -0.384
    # v0 = 0.073
    spot *= exp((r - q) * τ)
    df = exp(-r * τ)
    params= HestonParams(v0, κ, θ, ρ, σ)
    m = 200
    l = 8
    refPricer = makeCosCharFuncPricer(
        Complex,
        Float64,
        Float64(MathConstants.pi),
        params,
        τ,
        2048*4,
        32,
    )
    pricer =
        makeCosCharFuncPricer(Complex, Float64, Float64(MathConstants.pi), params, τ, m, l)
    refPrices = Vector{Float64}(undef, 0)
    prices = Vector{Float64}(undef, 0)
    for strike = 0.4:0.012:1.6
        refPrice = priceEuropean(refPricer, false, strike, spot, df)
        price = priceEuropean(pricer, false, strike, spot, df)
        println(
            strike,
            " ",
            price,
            " ",
            refPrice,
            " ",
            price - refPrice,
            " ",
            price / refPrice - 1,
        )
        push!(refPrices, refPrice)
        push!(prices, price)
    end
    rmse = rmsd(prices, refPrices)
    println("RMSE ", rmse)
    @test isless(rmse, 1e-4)
    @test isless(1e-5, rmse)
    #a somewhat large RMSE is expected, increase m for a workaround
    cinf = (params.v0+params.κ*params.θ*τ)/params.σ *sqrt(1-params.ρ^2)
    #exp(-cinf*umax) = eps
    #2*exp(-cinf*mmax*pi/(b-a))/(mmax*pi) = eps => exp(cinf*mmax*pi/(b-a))*mmax*pi/(b-a)*cinf = 2/eps/(b-a)*cinf => cinf*mmax*pi/(b-a) = lambertW( 2/eps/(b-a)*cinf)
    eps = 1e-8
    umax = -log(eps)/cinf
    mmax = ceil(Int,umax*(pricer.b-pricer.a)/pi)+1
    mmax = ceil(Int, lambertW(2*cinf/(eps*(pricer.b-pricer.a)))/(cinf*pi)*(pricer.b-pricer.a))
    pricer =
        makeCosCharFuncPricer(Complex, Float64, Float64(MathConstants.pi), params, τ, mmax, l)
        for (i, strike) = enumerate(0.4:0.012:1.6)
            refPrice = refPrices[i]
            price = priceEuropean(pricer, false, strike, spot, df)
            # println(
            #     strike,
            #     " ",
            #     price,
            #     " ",
            #     refPrice,
            #     " ",
            #     price - refPrice,
            #     " ",
            #     price / refPrice - 1,
            # )
            prices[i] = price
        end
        rmse = rmsd(prices, refPrices)
        println("RMSE ", rmse)
 end
#
# @testset "LFShortA" begin
#     r = 0.0
#     q = 0.0
#     κ = 0.254
#     θ = 0.320
#     σ = 0.344
#     ρ = -0.557
#     v0 = 0.826
#     τ = 0.0182
#     spot = 1000.0
#
#     spot *= exp((r - q) * τ)
#     df = exp(-r * τ)
#     params = HestonParams(v0, κ, θ, ρ, σ)
#     m = 200
#     l = 8
#     refPricer = makeCosCharFuncPricer(
#         Complex,
#         Float64,
#         Float64(MathConstants.pi),
#         params,
#         τ,
#         2048*4,
#         32,
#     )
#     pricer =
#         makeCosCharFuncPricer(Complex, Float64, Float64(MathConstants.pi), params, τ, m, l)
#     strike = 1400.0
#         refPrice = priceEuropean(refPricer, true, strike, spot, df)
#         price = priceEuropean(pricer, true, strike, spot, df)
#         println(
#             strike,
#             " ",
#             price,
#             " ",
#             refPrice,
#             " ",
#             price - refPrice,
#             " ",
#             price / refPrice - 1,
#         )
#
#     cinf = (params.v0+params.κ*params.θ*τ)/params.σ *sqrt(1-params.ρ^2)
#     #exp(-cinf*umax) = eps
#     #2*exp(-cinf*mmax*pi/(b-a))/(mmax*pi) = eps => exp(cinf*mmax*pi/(b-a))*mmax*pi/(b-a)*cinf = 2/eps/(b-a)*cinf => cinf*mmax*pi/(b-a) = lambertW( 2/eps/(b-a)*cinf)
#     eps = 1e-8
#     umax = -log(eps)/cinf
#     mmax = ceil(Int,umax*(pricer.b-pricer.a)/pi)+1
#     mmax = ceil(Int, sqrt(   lambertW(2*v0*τ/(eps*eps*(pricer.b-pricer.a)))/(v0*τ*pi)*(pricer.b-pricer.a)) )
#     mmax = max(32, ceil(Int, lambertW(2*cinf/(eps*(pricer.b-pricer.a)))/(cinf*pi)*(pricer.b-pricer.a)))
#     pricer =
#         makeCosCharFuncPricer(Complex, Float64, Float64(MathConstants.pi), params, τ, mmax, l)
#             price = priceEuropean(pricer, true, strike, spot, df)
#             println(
#                 strike,
#                 " ",
#                 price,
#                 " ",
#                 refPrice,
#                 " ",
#                 price - refPrice,
#                 " ",
#                 price / refPrice - 1,
#             )
# end
#
#
# @testset "LF1Y" begin
#     r = 0.0
#     q = 0.0
#     κ = 0.1
#     θ = 0.01
#     σ = 2.0
#     ρ = 0.5
#     v0 = 0.0225
#     τ = 1.0
#     spot = 1000000.0
#
#     spot *= exp((r - q) * τ)
#     df = exp(-r * τ)
#     params = HestonParams(v0, κ, θ, ρ, σ)
#     m = 200
#     l = 10
#     refPricer = makeCosCharFuncPricer(
#         Complex,
#         Float64,
#         Float64(MathConstants.pi),
#         params,
#         τ,
#         2048*4,
#         32,
#     )
#     pricer =
#         makeCosCharFuncPricer(Complex, Float64, Float64(MathConstants.pi), params, τ, m, l)
#     strike = spot*0.25
#         refPrice = priceEuropean(refPricer, false, strike, spot, df)
#         price = priceEuropean(pricer, false, strike, spot, df)
#         println(
#             strike,
#             " ",
#             price,
#             " ",
#             refPrice,
#             " ",
#             price - refPrice,
#             " ",
#             price / refPrice - 1,
#         )
#
#     cinf = (params.v0+params.κ*params.θ*τ)/params.σ *sqrt(1-params.ρ^2)
#     #exp(-cinf*umax) = eps
#     #2*exp(-cinf*mmax*pi/(b-a))/(mmax*pi) = eps => exp(cinf*mmax*pi/(b-a))*mmax*pi/(b-a)*cinf = 2/eps/(b-a)*cinf => cinf*mmax*pi/(b-a) = lambertW( 2/eps/(b-a)*cinf)
#     eps = 1e-10
#     umax = -log(eps)/cinf
#     mmax = ceil(Int,umax*(pricer.b-pricer.a)/pi)+1
#     mmax = ceil(Int, sqrt(   lambertW(2*v0*τ/(eps*eps*(pricer.b-pricer.a)))/(v0*τ*pi)*(pricer.b-pricer.a)) )
#     mmax = max(32, ceil(Int, lambertW(2*cinf/(eps*(pricer.b-pricer.a)))/(cinf*pi)*(pricer.b-pricer.a)))
#     pricer =
#         makeCosCharFuncPricer(Complex, Float64, Float64(MathConstants.pi), params, τ, mmax, l)
#             price = priceEuropean(pricer, false, strike, spot, df)
#             println(
#                 strike,
#                 " ",
#                 price,
#                 " ",
#                 refPrice,
#                 " ",
#                 price - refPrice,
#                 " ",
#                 price / refPrice - 1,
#             )
# end
@testset "LFShort" begin
    r = 0.0
    q = 0.0
    κ = 1.0
    θ = 0.1
    σ = 1.0
    ρ = -0.9
    v0 = 0.1
    τ = 2.0 / 365
    spot = 1.0

    spot *= exp((r - q) * τ)
    df = exp(-r * τ)
    params = HestonParams(v0, κ, θ, ρ, σ)
    m = 256
    l = 12
    refPricer = makeCosCharFuncPricer(
        Complex,
        Float64,
        Float64(MathConstants.pi),
        params,
        τ,
        1024,
        32,
    )
    pricer =
        makeCosCharFuncPricer(Complex, Float64, Float64(MathConstants.pi), params, τ, m, l)
    for strike = 1.0:0.025:1.5
        refPrice = priceEuropean(refPricer, true, strike, spot, df)
        price = priceEuropean(pricer, true, strike, spot, df)
        println(
            strike,
            " ",
            price,
            " ",
            refPrice,
            " ",
            price - refPrice,
            " ",
            price / refPrice - 1,
        )
        @test isapprox(refPrice, price, atol = 1e-15)
    end
end


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
    pricer =
        makeCosCharFuncPricer(Complex, Float64, Float64(MathConstants.pi), params, τ, m, l)
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
    pricer =
        makeCosCharFuncPricer(Complex, Float64, Float64(MathConstants.pi), params, τ, m, l)
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
