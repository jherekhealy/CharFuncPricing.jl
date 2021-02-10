using Nemo, CharFuncPricing, Test

#include("hestonfloat.jl")
include("schobelzhu.jl")
#Testing https://financepress.com/2019/02/15/heston-model-reference-prices/
@testset "HestonALSmallPrices" begin
    R = ArbField(512)
    CC = AcbField(512)
    r = R(QQ(0, 100))
    q = R(QQ(0, 100))
    κ = R(2)
    θ = R(QQ(1, 10000))
    σ = R(QQ(1, 10000))
    v0 = R(QQ(1, 10000))
    ρ = R(-0.5)
    τ = R(0.0025)
    strike = R(100)
    spots = [R("100.0001"), R(101), R(110), R(200), R(1000), R(10000)]
    params = HestonParams{arb}(v0, κ, θ, ρ, σ)

    l = 256 #48
    m = 1024 * 4 * 16
    cf = DefaultCharFunc{HestonParams{arb},acb,AcbField}(params, CC)
    pricer = makeCosCharFuncPricer(cf, τ, m, l)
    for spot in spots
        price = priceEuropean(pricer, false, strike, spot, R(1))
        println(spot, " ", price)
    end

    κ = big(2.0)
    θ = big(1.0)/big(10000.0)
    σ =  big(1.0)/big(10000.0)
    v0 =  big(1.0)/big(10000.0)
    ρ = big(-0.5)
    τ = big(0.0025)
    strike = big(100.0)
    spots = [BigFloat("100.0001"), big(101.0), big(110.0), big(200.0), big(1000.0), big(10000.0)]
    params = HestonParams{BigFloat}(v0, κ, θ, ρ, σ)
    cf = DefaultCharFunc{HestonParams{BigFloat},Complex{BigFloat},Type{Complex{BigFloat}}}(params, Complex{BigFloat})
    pricer = makeCosCharFuncPricer(cf, τ, m, l)
    for spot in spots
        price = priceEuropean(pricer, false, strike, spot, big(1.0))
        println(spot, " ", price)
    end
end

@testset "AlanSet" begin
    R = ArbField(256)
    CC = AcbField(256)
    r = R(QQ(1, 100))
    q = R(QQ(2, 100))
    κ = R(4.0)
    θ = R(QQ(1, 4))
    σ = R(1.0)
    v0 = R(QQ(4, 100))
    ρ = R(-0.5)
    τ = R(1.0)
    spot = R(100)
    strikes = [R(80), R(90), R(100), R(110), R(120)]
    alanPuts = [
        R("7.958878113256768285213263077598987193482161301733"),
        R("12.017966707346304987709573290236471654992071308187"),
        R("17.055270961270109413522653999411000974895436309183"),
        R("23.017825898442800538908781834822560777763225722188"),
        R("29.811026202682471843340682293165857439167301370697"),
    ]
    alanCalls = [
        R("26.774758743998854221382195325726949201687074848341"),
        R("20.933349000596710388139445766564068085476194042256"),
        R("16.070154917028834278213466703938231827658768230714"),
        R("12.132211516709844867860534767549426052805766831181"),
        R("9.024913483457835636553375454092357136489051667150"),
    ]
    #refPuts = [R("7.95887811325676828521326060761429303089865693725960319205094095681790030"), R("12.01796670734630498770957286168505994488766901674359014792496983503358626"), R("17.05527096127010941352265395975862503115076404332454551322405076347166800"),R("23.01782589844280053890878183747069958885882262385534095209843808561556738"), R("29.81102620268247184334068229213413372849045238688265116546863785487435207")]
    #refCalls = [R("26.77475874399885422138219285574225503910357048386701361386103443917276629"),R("20.93334900059671038813944533801265637537179175081262590089627586445697497"), R("16.07015491702883427821346666428585588391409596485520659735656933996357944"), R("12.13221151670984486786053477019756486390136373284762736739216920917600154"), R("9.02491348345783563655337545306063342581220268333656291192358152550330896")]
    spot *= exp((r - q) * τ)
    df = exp(-r * τ)
    params = HestonParams{arb}(v0, κ, θ, ρ, σ)

    l = 48
    m = 1024 * 4
    cf = DefaultCharFunc{HestonParams{arb},acb,AcbField}(params, CC)
    pricer = makeCosCharFuncPricer(cf, τ, m, l)
    for (strike, refCall, refPut) in zip(strikes, alanCalls, alanPuts)
        price = priceEuropean(pricer, false, strike, spot, df)
        println(Float64(strike), " P ", price, " ", price - refPut)
        #push!(refPuts,price)
        @test isapprox(Float64(price - refPut), 0, atol = 1e-23)
        price = priceEuropean(pricer, true, strike, spot, df)
        println(Float64(strike), " C ", price, " ", price - refCall)
        @test isapprox(Float64(price - refCall), 0, atol = 1e-23)
        #push!(refCalls,price)
    end
end


@testset "AlanJoshiSet" begin
    R = ArbField(256)
    CC = AcbField(256)
    r = R(QQ(1, 100))
    q = R(QQ(2, 100))
    κ = R(4.0)
    θ = R(QQ(1, 4))
    σ = R(1.0)
    v0 = R(QQ(1, 100))
    ρ = R(-0.5)
    τ = R(QQ(1, 100))
    spot = R(100)
    strikes = [R(90), R(95), R(100), R(105), R(110)]
    alanPuts = [
        R("4.5183603586861772614990106188215872180542e-8"),
        R("0.000461954855653851579672612557018857858641926937"),
        R("0.477781171629504680023239655436072890669645669297"),
        R("5.009501052563650299130635110520904481889436667608"),
        R("10.008998550115123724684210555728039829315964456261"),
    ]
    alanCalls = [
        R("9.989001595065276544935948045293485530832966049263"),
        R("4.989963479738160122154264702582719627807098780529"),
        R("0.467782671512844263098248405184095087949465507760"),
        R("2.527447823194706060519991248106500619490942e-6"),
        R("1.29932760052624920704881258510264466e-13"),
    ]
    refPuts = [
        R("4.518360358686177267846283260505342170483575627777422720994e-8"),
        R("0.00046195485565385157967261251100128843593583885256695108090334577"),
        R("0.47778117162950468002323965543594159177844084992306552196870005758"),
        R("5.00950105256365029913063511052075051116901889480844354030472812057"),
        R("10.00899855011512372468421052192118981283250808930047516269268570714"),
    ]
    refCalls = [
        R("9.98900159506527654493594810876621194767051557355651881698849965467"),
        R("4.989963479738160122154264702536702058384392692444362488921902796131"),
        R("0.46778267151284426309824840518396378905826068838589405843642651357"),
        R("2.52744782319470606051999109413578020171814230507539918158220e-6"),
        R("1.2993276005259111385486477505389750536969641386617440e-13"),
    ]
    spot *= exp((r - q) * τ)
    df = exp(-r * τ)
    params = HestonParams{arb}(v0, κ, θ, ρ, σ)
    l = 48
    m = 1024 * 4
    pricer = makeCosCharFuncPricer(
        DefaultCharFunc{HestonParams{arb},acb,AcbField}(params, CC),
        τ,
        m,
        l,
    )
    for (strike, refCall, refPut) in zip(strikes, alanCalls, alanPuts)
        price = priceEuropean(pricer, false, strike, spot, df)
        println(Float64(strike), " P ", price, " ", price - refPut)
        #push!(refPuts,price)
        @test isapprox(Float64(price - refPut), 0, atol = 1e-25)
        price = priceEuropean(pricer, true, strike, spot, df)
        println(Float64(strike), " C ", price, " ", price - refCall)
        @test isapprox(Float64(price - refCall), 0, atol = 1e-25)
        #push!(refCalls,price)
    end
end

@testset "ZeroKappa" begin
    R = ArbField(256)
    CC = AcbField(256)
    κ = R(0.0)
    θ = R(QQ(1, 4))
    σ = R(1.0)
    v0 = R(QQ(4, 100))
    ρ = R(-0.5)
    τ = R(1.0)
    params = HestonParams{arb}(v0, κ, θ, ρ, σ)
    cumulantsZ = computeCumulants(params, τ)
    κ = R(QQ(1, 100))
    params = HestonParams{arb}(v0, κ, θ, ρ, σ)
    cumulants = computeCumulants(params, τ)

    for (cz, c) in zip(cumulantsZ, cumulants)
        println(cz, " ", c)
        @test isapprox(Float64(cz), Float64(c), rtol = 1e-1)
    end
end
