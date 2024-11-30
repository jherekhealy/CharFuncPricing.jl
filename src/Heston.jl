export HestonParams, evaluateCharFunc, evaluateLogCharFunc, computeCumulants
export BlackParams, evaluateCharFuncAndDerivative
export HestonCVCharFunc, cinf

struct HestonParams{T}
    v0::T
    κ::T
    θ::T
    ρ::T
    σ::T
end
DefaultCharFunc(params::HestonParams{Float64}) =
    DefaultCharFunc{HestonParams{Float64},Complex{Float64}}(params)

DefaultCharFunc(params::HestonParams{T}) where {T} =
    DefaultCharFunc{HestonParams{T},Complex{T}}(params)

HestonCVCharFunc(heston::CharFunc{HestonParams{T},CR}, τ, kind::CVKind=InitialControlVariance()) where {T,CR} =
    CVCharFunc{HestonParams{T},BlackParams{T},CR}(
        heston,
        DefaultCharFunc{BlackParams{T},CR}(
            BlackParams{T}(computeControlVariance(heston, τ, kind))
        )
    )

makeCVCharFunc(heston::CharFunc{HestonParams{T},CR}, τ, kind::CVKind) where {T,CR} = HestonCVCharFunc(heston, kind)


@inline function computeControlVariance(
    cf::CharFunc{HestonParams{TT}},
    τ::T, kind::FullControlVariance
)::T where {T,TT}
    p = model(cf)
    ektm = 1 - exp(-p.κ * τ)
    return ektm * (p.v0 - p.θ) / (p.κ * τ) + p.θ
end

@inline function computeControlVariance(
    cf::CharFunc{HestonParams{TT}},
    τ::T, kind::InitialControlVariance
)::T where {T,TT}
    p = model(cf)
    p.v0
end


#cinf(params::HestonCVParams{T}, τ::T) where {T} = cinf(params.heston, τ)

cinf(params::HestonParams{T}, τ::T) where {T} =
    (params.v0 + params.κ * params.θ * τ) / params.σ * sqrt(1 - params.ρ^2)

#CC the Nemo arbField or Complex type. z: Nemo complex number or Complex number.
@inline function evaluateCharFunc(
    p::CharFunc{MT,CR},
    z::CT,
    τ::T,
) where {T,CR,CT,MT}
    return exp(evaluateLogCharFunc(p, z, τ))
end

@inline evaluateLogCharFunc(
    p::CharFunc{HestonParams{T},CR},
    z::CT,
    τ::T,
) where {T,CR,CT} = evaluateLogCharFuncAL(p, z, τ)


@inline function evaluateLogCharFuncAL(
    cf::CharFunc{HestonParams{T},CR},
    z::CT,
    τ::T,
) where {T,CR,CT}
    #follows Andersen-Lake fast implementation
    p = model(cf)
    v0 = p.v0
    κ = p.κ
    θ = p.θ
    ρ = p.ρ
    σ = p.σ
    cc1 = oneim(cf)
    α = -(z * (z + cc1)) * σ^2
    β = κ - σ * ρ * z * cc1
    D = sqrt(β^2 - α)
    rr = (real(β) * real(D) + imag(β) * imag(D)) > 0 ? α / (β + D) : (β - D)
    y = D == 0 ? -τ / 2 : expm1(-D * τ) / (2 * D)
    l = log1p(-rr * y)
    A = (rr * τ - 2 * l) * (κ * θ / σ^2)
    B = z * (z + cc1) * y / (1 - rr * y)
    return A + B * v0
end

#Faster evaluation of the characteristic function phi(ix) where x is real. Useful for optimal alpha.
@inline function evaluateLogCharFuncAtImag(
    cf::CharFunc{HestonParams{T},CR},
    imz::T,
    τ::T,
)::T where {T,CR}
    #follows Andersen-Lake fast implementation
    p = model(cf)
    v0 = p.v0
    κ = p.κ
    θ = p.θ
    ρ = p.ρ
    σ = p.σ
    cc1 = oneim(cf)
    α = (imz * (imz + one(T))) * σ^2
    β = κ + σ * ρ * imz
    D2 = β^2 - α
    local ch, sh
    if D2 >= 0
        D = sqrt(D2)
        ch = cosh(D * τ / 2)
        sh = (D == 0) ? τ / 2 : sinh(D * τ / 2) / D
    else
        D = sqrt(-D2)
        sh, ch = sincos(D * τ / 2)
        sh /= D
    end
    A = κ * θ / σ^2 * (β * τ - log((ch + β * sh)^2))
    B = imz * (imz + 1) * sh / (ch + β * sh)
    return A + B * v0
end

@inline function evaluateCharFuncAndDerivative(
    p::CharFunc{MAINT,CR},
    z::CT,
    τ::T,
)::Tuple{CR,CR} where {T,CR,CT,MAINT}
    arg, arg_d = evaluateLogCharFuncAndDerivative(p, z, τ)
    phi = exp(arg)
    phi_d = arg_d * phi
    return phi, phi_d
end

@inline function evaluateLogCharFuncAndDerivative(
    cf::CharFunc{HestonParams{T},CR},
    z::CT,
    τ::T,
)::Tuple{CR,CR} where {T,CR,CT}
    #derivative towards real part of z
    p = model(cf)
    v0 = p.v0
    κ = p.κ
    θ = p.θ
    ρ = p.ρ
    σ = p.σ
    cc1 = oneim(cf)
    α = -(z * (z + cc1)) * σ^2
    α_d = -(2 * z + cc1) * σ^2
    β_d = -σ * ρ * cc1
    β = κ + β_d * z

    D = sqrt(β^2 - α)
    D_d = if D == zero(T)
        zero(T)
    else
        (2 * β * β_d - α_d) / (2 * D)
    end
    local rr
    local rr_d
    if (real(β) * real(D) + imag(β) * imag(D)) > 0  
        rr = α / (β + D)
        rr_d = α_d / (β + D) - rr * (β_d + D_d) / (β + D)
    else
        rr = (β - D)
        rr_d = β_d - D_d
    end
    local y
    local y_d
    if D != 0
        em = expm1(-D * τ)
        y = em / (2 * D)
        y_d = (y * (-2 * D_d * (D * τ + 1)) + (-D_d * τ)) / (2 * D) # ((em+1)*(-D_d * τ)*(2*D) - em*(2*D_d))/(2*D)^2
    else
        y = -τ / 2
        y_d = zero(T)
    end
    l = log1p(-rr * y)
    rry_d = (rr_d * y + rr * y_d)
    l_d = -rry_d / (1 - rr * y)
    Bnum = z * (z + cc1) * y
    Bdenom = (1 - rr * y)
    B = Bnum / Bdenom
    B_d = ((2 * z + cc1) * y + y_d * z * (z + cc1)) / Bdenom + rry_d * B / Bdenom
    #if rr*y == one(T)       
    #end
    A = (rr * τ - 2 * l) * (κ * θ / σ^2)
    A_d = (rr_d * τ - 2 * l_d) * (κ * θ / σ^2)
    return A + B * v0, A_d + B_d * v0
end

struct CuiCharFunc{MT,CR} <: CharFunc{MT,CR} #model type, return type (e.g. Complex or acb),
    model::MT
end

CuiCharFunc(params::HestonParams{Float64}) =
    CuiCharFunc{HestonParams{Float64},Complex,Type}(params, Complex)
@inline model(cf::CuiCharFunc) = cf.model
@inline field(cf::CuiCharFunc) = cf.field
@inline evaluateLogCharFunc(
    p::CuiCharFunc{HestonParams{T},CR},
    z::CT,
    τ::T,
) where {T,CR,CF,CT} = evaluateLogCharFuncAL(p, z, τ)


function evaluateLogCharFuncCui(
    cf::CharFunc{HestonParams{T},CR},
    z::CT,
    τ::T,
)::CR where {T,CR,CT}
    p = model(cf)

    v0 = p.v0
    κ = p.κ
    θ = p.θ
    ρ = p.ρ
    σ = p.σ
    cc1 = oneim(cf)
    iu = cc1 * z
    ξ = κ - σ * ρ * iu
    α = iu * (-iu + 1)
    d = sqrt(ξ^2 + σ^2 * α)
    dht = d * (τ / 2)
    sh1 = sinh(dht)
    ch1 = cosh(dht)
    A2v = (d * ch1 + ξ * sh1)
    A1 = α * sh1
    # edt = ch1 - sh1
    # logB = (κ * τ / 2 - dht) - log1p( (ξ - d) * sh1/ d * edt )

    edt = (ch1 + sh1)
    logB = (κ * τ / 2 - dht) - log(A2v / (d * edt))
    return -κ * θ * ρ * τ / σ * iu - A1 / A2v * v0 + 2 * κ * θ / σ^2 * logB
end

function evaluateFourthDerivative(cf::CharFunc{HestonParams{T},CR}, τ::T, u::T) where {T,CR}
    m = model(cf)
    e1 = (0 + 1im) * u
    e3 = e1 * m.σ * m.ρ
    e4 = e3 - m.κ
    e5 = m.σ^2
    e7 = e4^2 - u * (0 - 1im - u) * e5
    e8 = sqrt(e7)
    e13 = 2 * ((0 + 1im) * e4 * m.ρ) - (0 - 1im - 2 * u) *
                                       m.σ
    e14 = (0 + 1im) * m.ρ
    e16 = 0.5 * (e13 / e8)
    e18 = m.κ - (e3 + e8)
    e20 = m.κ + e8 - e3
    e21 = e16 - e14
    e23 = exp(-(τ * e8))
    e24 = m.ρ^2
    e26 = 2 * ((-1 + 0im) * e24)
    e27 = 2 + e26
    e28 = e13^2
    e29 = e28 / e7
    e30 = 0.5 * e29
    e31 = e27 - e30
    e32 = e14 + e16
    e33 = 1 - e23 * e18 / e20
    e34 = τ * e18
    e35 = 0.5 * e34
    e36 = 0.5 + e35
    e37 = τ * e13
    e38 = e31 / e8
    e39 = e37 / e8
    e40 = 0.5 * e39
    e41 = 0.5 * e38
    e42 = e21 / e20
    e43 = e42 + e40
    e45 = e21 * e18 / e20
    e46 = e36 * e13
    e48 = e21^2 / e20
    e51 = e43 * e18 + e14 + e16
    e52 = e41 - e48
    e54 = e45 + e46 / e8 + e14
    e55 = 1 - e23
    e56 = e32 * e21
    e57 = e52 * e18
    e58 = τ * e32
    e59 = (e33 * e5)^2
    e61 = 1 - e18 / e20
    e63 = 2 * e27 - e29
    e65 = e45 + e14 + e16
    e66 = e46 / e7
    e67 = 0.5 * e66
    e68 = τ * e31
    e69 = e33 * e20
    e70 = 0.5 * e31
    e71 = e57 - e56
    e72 = e68 / e8
    e73 = e58 * e13
    e74 = m.σ^4
    e75 = 0.5 * e63
    e76 = 0.5 * e73
    e77 = e71 / e20
    e78 = 0.5 * e72
    e79 = e75 + e70
    e80 = e36 * e31
    e81 = e7 * e8
    e83 = e52 / e20 + e78
    e84 = 0.5 * e54
    e85 = e51 * e23
    e86 = e36 * e27
    e87 = e43 * e32
    e88 = e83 * e18
    e90 = e51 * e33 * e23
    e91 = 0.5 * e55
    e92 = e90 * e74
    e94 = e59 * e21 + 2 * e92
    e95 = e79 * e13
    e96 = e84 + 0.5 * e32
    e97 = e80 - e76
    e98 = τ / e8
    e99 = e32 * e13
    e100 = 0.5 / e7
    e101 = e32 * e55
    e103 = 0.5 * e23 - e91
    e104 = τ * e96
    e105 = e95 / e81
    e106 = e97 / e8
    e107 = e51 * e54
    e108 = e88 + e41
    e109 = e65 / e61
    e110 = 0.5 * e105
    e111 = e37 * e23
    e113 = e107 * e23 / e33
    e114 = e13 * e31
    e116 = 0.5 * e98 + e100
    e117 = e65^2
    e118 = e117 / e61
    e119 = e77 + e106
    e120 = e38 - e48
    e122 = e86 - (e67 + 0.5 * e58) * e13
    e124 = 0.5 * (e111 * e18 / e8) - e101
    e127 = e120 * e21 / e20
    e128 = e57 - (e113 + e56)
    e129 = τ * e54
    e130 = e127 + e110
    e131 = e116 * e28
    e134 = e54 * e23 / e33 - e109
    e135 = e130 * e18
    e137 = e32 * e52
    e138 = e57 - (e118 + e56)
    e139 = e114 / e7
    e140 = e94 * e54
    e141 = e33 * m.σ
    e144 = e40 - e51 * e55 / e69
    e145 = e122 / e8
    e146 = e59 * e20
    e147 = e128 / e20
    e148 = 0.5 * e139
    e149 = e27 - e131
    e150 = e108 - e87
    e152 = (-4 + 0im) * e24 + 4
    e157 = 0.5 * (e149 * e18 - e99) - 0.5 * e99
    e159 = 0.5 * (e55 * e31)
    e160 = e140 / e59
    e162 = e138 / e20 + e41
    e163 = τ * e157
    e164 = e162 / e61
    e166 = e141^2
    e167 = e57 - (e160 + e56)
    e170 = e129 * e103 * e13 / e8
    e171 = e27 - e29
    e173 = e54 * e32 * e55
    e174 = e167 / e20
    e176 = e21 * e33
    e177 = e36 * e171
    e178 = e86 - (e67 + e104) * e13
    e179 = e152 - e30
    e181 = e163 * e23 - e159
    e182 = e77 + e145
    e183 = e177 - e76
    e184 = e80 - e104 * e13
    e186 = e134 * e21 / e20
    e188 = 0.0
    e189 = e51 * e124
    e190 = e184 / e8
    e191 = 0.5 * e183
    e192 = e179 * e36
    e193 = e51 * e21
    e194 = e85 / e20
    e195 = e67 + e58
    e196 = τ * e51
    e197 = (0 + 1im) * e188
    e198 = e196 * e13
    e199 = e186 + e164
    e203 = e51^2 * e23 / e69 + e87
    e206 = e144 * e23 * e18 - e101
    e207 = 0.5 * (τ * e33 * e13 / e8)
    e212 = e85 + e176
    e214 = e85 / e33 - e109
    e216 = 0.5 * (e114 / e8)
    e218 = 2 * ((e194 - e207) * e51 + e150 * e33) + 2 *
                                                    (e193 * e33 / e20)
    e219 = e147 + e190
    e226 = e218 * e23 * e74 + 0.5 * (e59 * e31 / e8)
    e227 = 0.25 * e38
    e228 = 0.5 * (e32 * e27 + e216)
    e229 = e192 - (e67 + 2 * e104) * e13
    e232 = e135 + e71 * e21 / e20 + e137 + e148
    e233 = e147 + e178 / e8
    e234 = e108 - e203
    e235 = 0.5 * e119
    e236 = e232 / e20
    e238 = e219 * e23 / e33
    e239 = e182 * e55
    e240 = e234 * e54
    e241 = e195 * e27
    e242 = 0.5 * e122
    e243 = e229 / e8
    e245 = e212 * e51 / e69
    e246 = e238 - e199
    e247 = e119 * e103
    e248 = e54 * e226
    e249 = e174 + e145
    e250 = e77 + e41
    e251 = (e191 + e242) / e7
    e252 = e181 / e8
    e253 = e129 * e13
    e255 = τ * e23 / e8
    e256 = e189 / e166
    e257 = e85 / e69
    e259 = e128 * e21 / e20
    e260 = 0.25 * e72
    e262 = 0.25 * e98 + e100
    e263 = 0.5 * e178
    e264 = 2 * e250
    e265 = τ * e79
    e266 = e245 + e87
    e268 = e257 + e40
    e270 = e249 * e55 + e170
    e272 = e119 * e32 + 0.5 * (e54 * e31 / e8)
    e274 = (e251 + e260) * e13 + e241
    e275 = e54 * e103
    e283 = (e57 - (e56 + 2 * e118)) / e20 + e41 + e264
    e286 = (e124 / e141 - e54 * e55 * e23 * e18 * m.σ^3 / e146) *
           m.v0 + e197 - (2 * (e134 / e20) + e58) * m.κ *
                         m.θ / m.σ
    e290 = e103 / e7
    e292 = e14 + e84 + e16
    e293 = 0.5 * e255
    e294 = 2 * e173
    e296 = 2 * (e116 * e27) - e262 * e28 / e7
    e298 = e198 * e23 / e8
    e300 = e265 * e13 / e81
    e301 = e152 - e131
    e303 = e270 * e18 - e173
    e304 = e272 * e55
    e305 = e274 / e8
    e306 = e150 * e124
    e308 = e138 * e21 / e20
    e309 = e166 * e20
    e312 = e296 * e18 + e41
    e314 = 0.5 * e290
    e315 = e70 - 0.5 * (e144 * e13)
    e316 = 2 * e56
    e317 = 2 * e57
    e318 = e111 / e8
    e321 = ((e283 * e65 / e61 + e135 + e308 + e137 + e148) / e20 +
            e110) / e61
    e326 = e51 * e206 / e69
    e328 = e119 * e94 + e248
    e329 = e239 + e170
    e331 = (e206 * m.v0 / e33 - (2 * (e214 / e20) + e58) *
                                m.κ * m.θ) / m.σ + e197
    e334 = e32 * e144
    e335 = (e67 + τ * e292) * e27
    e336 = e176 + 2 * e85
    e337 = (e70 - 0.5 * e198) / e8
    e341 = e227 + e235
    e343 = 0.5 * (e301 * e32 + e312 * e13) + e228
    e344 = 0.5 * (e134 * e31 / e8)
    e345 = e314 + e293
    e346 = 0.5 * (e36 * e63)
    e347 = e191 + e263
    e349 = e159 + 0.5 * (e73 * e23)
    e350 = 0.5 * e298
    e351 = 0.5 * e300
    e353 = 0.5 * (e68 * e23 / e8) - 0.5 * (e55 * e63 / e7)
    e354 = 2 * (e94 * e55 / e146)
    e355 = e23 - e91
    e357 = τ * e315 / e8
    e360 = e303 * e5 / e59 + e256
    e362 = e214 * e21 / e20
    e364 = e328 * e55 + e140 * (0.5 * e318 - e354)
    e366 = e329 * e18 - e294
    e367 = e236 + e305
    e368 = e326 + e334
    e370 = e189 * e23 / e69
    e371 = e182 * e23
    e372 = e54 * e345
    e373 = e174 + e106
    e374 = e88 + e337
    e376 = e347 / e7
    e377 = e349 / e8
    e378 = e33 * e8
    e379 = e343 + 0.5 * (e163 * e13 / e8)
    e380 = 0.5 * e233
    e381 = 0.5 * (e32 * e31)
    e382 = e192 - e195 * e13
    e383 = e243 + (e317 - (e113 + e316)) / e20
    e384 = e243 + 2 * e147
    e385 = (e364 * e18 + e366 * e94) / e146
    e386 = e367 * e55
    e392 = (e51 * (e77 + e190) + e240) * e23 / e33 + e135 +
           e259 + e137 + e148
    e393 = e362 + e164
    e394 = e373 * e55
    e395 = e372 * e13
    e396 = e275 * e27
    e397 = e33 * e13
    e399 = e227 + e380 + e235
    e400 = 0.5 * e371
    e401 = 2 * e164
    e403 = e253 * e355 / e8
    e405 = τ * e379 * e23
    e407 = e382 / e8 + 2 * e174
    e409 = e383 * e51 + e240
    e412 = e384 * e23 / e33 - (2 * e186 + e401)
    e413 = e392 / e20
    e415 = (e394 + e170) * e18 - e173
    e417 = (e247 + e400 - e395) * e13 + e396
    e418 = ((e21 * (1.5 * e38 - 2 * e48) / e20 + e110) / e20 +
            e351) * e18
    e421 = e233 * e23 / e33 - e199
    e422 = e239 + e403
    e424 = e309^2
    e427 = ((e346 + 0.5 * e97) * e13 / e7 + τ * (e381 +
                                                 e228)) / e8
    e430 = e108 - e266
    e431 = (0.5 * (e43 * e31) + 0.5 * (e95 / e7)) / e8
    e433 = e79 * e27 - ((0.5 * (e7 / e8) + e8) * e79 / e8 +
                        0.75 * e63) * e28 / e7
    e434 = e27 * e31
    e436 = e75 + 2 + e26
    e437 = 2 * e87
    e438 = 2 * (e83 * e32)
    e439 = 2 * e88
    e440 = 0.0
    e442 = e407 * e55 + 2 * e170
    e448 = e409 * e23 / e33 + e135 + e259 + e137 + e148
    e450 = e412 * e21 + e344
    e452 = (((e266 - e108) * e55 - e350) / e69 + e357) *
           e18 - e368
    e458 = e268 * e51
    e459 = e246 * e21
    e461 = (((e57 - (e56 + 2 * (e107 * e33 * e23 * e74 / e59))) / e20 +
             e106) * e94 + e248) / e59
    e462 = (e430 * e55 + e350) / e69
    e463 = e331 * e286
    e464 = e236 + e427
    e465 = e51 * e336
    e466 = e193 / e20
    e469 = e51 * e181 / e8 + e306
    e472 = e85 * e8 / e20 + 0.5 * (e397 / e8)
    e474 = ((e252 - e370) / e33 - e415 * e23 * e74 / e146) *
           m.v0 - (e78 + 2 * (e246 / e20)) * m.κ * m.θ
    e476 = (e376 + τ * e399) * e13 + e335
    e477 = e275 * e31
    e482 = e38 + e439
    e484 = (τ * e417 / e8 - e386) * e18 - (e385 + e304 +
                                           e422 * e32)
    e485 = 0.5 * (e433 / e81)
    e486 = 0.5 * (e353 * e13)
    e487 = 0.5 * (e253 / e8)
    e488 = 2 * e135
    e489 = 2 * e137
    e490 = exp((e55 * e18 * m.v0 / e33 + (e34 - 2 *
                                                (log(e33) - log(e61))) * m.κ * m.θ) / e5 +
               e1 * e188)
    e492 = e129 * e32 * e13
    e493 = τ * e341
    e495 = (-6 + 0im) * e24 + 6
    e497 = e442 * e18 - e294
    e499 = e448 / e20 + e476 / e8
    e500 = e450 / e20
    e503 = (e452 * e23 - e377) * m.v0 / e33 - (e78 +
                                               2 * (((e108 - (e458 + e87)) * e23 / e33 - e393) / e20)) *
                                              m.κ * m.θ
    e504 = (e461 + e135 + e167 * e21 / e20 + e137 + e148) / e20
    e507 = e469 / e166 + e484 * e5 / e59 - 0.5 * (τ * e360 *
                                                  e13 / e8)
    e508 = e459 + e344
    e510 = ((e376 + e493) * e13 + e335) / e8
    e512 = e418 + e431 + e438
    e516 = e465 * e124 * e33 * e5 / e424
    e524 = (((e357 - e462) * e18 - e368) * e23 - e377) *
           m.v0 / e33 - (e78 + 2 * (((e374 - e203) * e23 / e33 -
                                     e393) / e20)) * m.κ * m.θ
    e525 = ((0.5 * (e120 * e31 / e8) - (e105 + 2 * e127) *
                                       e21) / e20 + e485) * e18
    e529 = 0.5 * e182
    e530 = 0.5 * ((e434 - (e436 - e30) * e28 / e7) / e7)
    e531 = 0.5 * (e52 * e31 / e8)
    e532 = 0.5 * e353
    e533 = 0.5 * (τ * e179 * e32)
    e534 = 2 * (e189 * e33 * e23 * e5 / e309)
    e535 = 2 * e246
    e536 = 2 * e421
    e537 = 2 * (e130 * e32)
    e538 = 2 * e252
    e539 = 3 * e252
    e541 = (e495 - e29) * e36
    e543 = e500 + e499 * e23 / e33 - e321
    e546 = e497 * e5 / e59 + e256
    e547 = e503 * e286
    e549 = (e413 + e268 * e219 + ((e346 + 0.5 * e184) *
                                  e13 / e7 + τ * (e341 * e13 + e96 * e27 + e381)) / e8) *
           e23 / e33
    e551 = (e413 + e268 * e233 + e510) * e23 / e33
    e553 = e507 / e20 - e516
    e554 = (e472 * e181 / e378 + e486 + e405) / e378
    e555 = e508 / e20
    e560 = (e51 * (e252 - e534) + e306) / e166
    e562 = (e194 - 0.5 * (e397 / e7)) * e13 + e33 * e27
    e564 = ((((e436 - e29) * e36 + 0.5 * e229 + e263 -
              e76) / e7 + τ * (e399 + 2 * e341)) * e13 + (e66 +
                                                          τ * (e292 + 2 * e96)) * e27 + e533) / e8
    e565 = e524 * e286
    e566 = (e336 / e69 + e40) * e124
    e567 = e134 * e13
    e569 = e51 * e13
    e570 = e150 * e13
    e572 = e54 * e13
    e573 = e374 - (e466 + e87)
    e574 = ((0.5 * (e181 / e7) + e532) * e13 + e405) / e8
    e577 = e59 * e13
    e578 = (e55 * (2 * e245 + e437 - e482) - e298) / e69
    e579 = (0.5 * (e492 / e8) + 2 * (e303 * e51 * e33 *
                                     e74 / e146)) * e23
    e586 = 0.5 * e249
    e587 = 0.5 * (e183 / e7)
    e588 = e486 + e405
    e589 = 0.5 * (e46 * e63 / e7 + τ * (0.5 * (e32 * e171) +
                                        e228))
    e590 = 0.5 * (e33 * e31 / e8)
    e591 = 0.5 * (e13 / e7)
    e592 = 1.5 * e58
    e593 = 1.5 * e72
    e594 = 2 * (e304 + 0.5 * (e492 * e23 / e8))
    e595 = 2 * e113
    e596 = e535 + e536
    e597 = 2 * e247
    e598 = e317 - (e160 + e316)
    e599 = 3 * e173
    e600 = 3 * e118
    e601 = 3 * e174
    e602 = 3 * e56
    e603 = 3 * e57
    e604 = 3 * e170
    ((((((((((((e466 + e88 + e337 - e87) * e23 + e590) *
              e51 + e212 * e150) * e55 + (e51 * (e318 - 2 * (e212 *
                                                             e55 / e69)) + e150 * e55) * e212) / e69 + e512 * e55 -
            τ * (0.5 * (e51 * e149 + e570) + 0.5 * e570) *
            e23 / e8) / e69 - τ * ((0.25 * e63 + 0.5 * e315) *
                                   e13 / e7 + 0.5 * (e144 * e27 + (e78 - e462) * e13)) / e8) *
          e18 - (((((e578 + 2 * e357) * e18 - (e326 + 2 *
                                                      e334)) * e23 - 2 * e377) * e51 + e430 * e206) / e69 +
                 (e578 + τ * (e27 - (e591 + 0.5 * e144) * e13) / e8) *
                 e32 + (0.5 * (e144 * e31) + 0.5 * (τ * e452 *
                                                    e13)) / e8)) * e23 - ((e532 - 0.5 * (e349 / e7)) *
                                                                          e13 + 0.5 * (τ * (e32 * (e27 - 0.5 * (τ * e28 / e8)) +
                                                                                            e216) * e23)) / e8) * m.v0 / e33 + (e351 + 2 *
                                                                                                                                       (((e268 * (e482 - (e458 + e437)) + ((e374 - e266) *
                                                                                                                                                                           e23 / e69 + e78) * e51 + e418 + e431 + e438) *
                                                                                                                                         e23 / e33 + (e214 * e52 + (((e27 - (e591 + 0.5 *
                                                                                                                                                                                    e196) * e13) / e8 + e439 - (e51 * (e40 + 2 *
                                                                                                                                                                                                                             e257) + e437)) * e23 / e33 - (e362 + e401)) *
                                                                                                                                                                   e21) / e20 - e321) / e20)) * m.κ * m.θ) *
       e286 + ((0.5 * (τ * e433 / e81) + 2 * (((((e383 *
                                                  e150 + e119 * e234 - ((((((((-10 + 0im) * e24 + 10 -
                                                                              1.5 * e29) * e36 - (e66 + 5 * e104) * e13) / e8 +
                                                                            (5 * e57 - (e595 + 5 * e56)) / e20) * e51 + 3 * e240) *
                                                                          e23 / e33 + e135 + e137 + e21 * (e603 - (e595 +
                                                                                                                   e602)) / e20 + 1.5 * e139 + 2 * (e135 + e137)) / e20 +
                                                                         e564) * e51 + (e418 + e51 * (2 * e150 - (e212 / e69 +
                                                                                                                  e40) * e51) * e23 / e69 + e431 + e438) * e54 + 0.5 *
                                                                                                                                                                 (τ * e409 * e13 / e8))) * e23 / e33 + (e128 * e52 -
                                                                                                                                                                                                        ((((e541 - (e67 + 3 * e104) * e13) / e8 + (e603 -
                                                                                                                                                                                                                                                   (e113 + e602)) / e20) * e51 + 2 * e240) * e23 / e33 +
                                                                                                                                                                                                         e259 + e139 + e488 + e489) * e21) / e20 + e525 +
                                                 e530 + e531 - e537) / e20 + (((e177 + e263 - e76) / e7 +
                                                                               τ * (e77 + ((e35 + 1.25) * e31 - e76) / e8 + e380)) *
                                                                              e27 - ((e347 * e13 / e7 + 0.5 * e476 + 0.5 * ((e587 +
                                                                                                                             e493) * e13 + e335) + e589) / e7 + τ * (0.25 * e105 +
                                                                                                                                                                     0.5 * e499 + 0.5 * (e413 + e510) + 0.5 * e464)) *
                                                                                    e13) / e8) * e23 / e33 + ((0.5 * (e412 * e31) + 0.5 *
                                                                                                                                    ((e238 - (e164 + 0.5 * (e567 / e7))) * e31 - 0.5 *
                                                                                                                                                                                 (e567 * e63 / e7))) / e8 - ((e384 * e268 + (e448 +
                                                                                                                                                                                                                             2 * e392) / e20 + e564) * e23 / e33 + (2 * e450 + 2 *
                                                                                                                                                                                                                                                                               e508) / e20 - 3 * e321) * e21) / e20 - (((e162 * e283 -
                                                                                                                                                                                                                                                                                                                         (((e283 + 2 * (e264 - e117 / (e61 * e20))) * e65 / e61 +
                                                                                                                                                                                                                                                                                                                           e21 * (e317 - (e316 + e600)) / e20 + e139 + e488 +
                                                                                                                                                                                                                                                                                                                           e489) / e20 + e105 + 2 * (e236 + e110)) * e65) / e61 +
                                                                                                                                                                                                                                                                                                                        (e138 * e52 - (((e57 - (e56 + e600)) / e20 + e41 +
                                                                                                                                                                                                                                                                                                                                        4 * e250) * e65 / e61 + e308 + e139 + e488 +
                                                                                                                                                                                                                                                                                                                                       e489) * e21) / e20 + e525 + e530 + e531 - e537) / e20 +
                                                                                                                                                                                                                                                                                                                       e485) / e61) / e20)) * m.κ * m.θ - (((((e573 *
                                                                                                                                                                                                                                                                                                                                                               e8 + 0.5 * (e569 / e8)) * e23 / e20 + (0.5 * e562 -
                                                                                                                                                                                                                                                                                                                                                                                                      2 * (e472^2 / e33)) / e8) * e181 - 2 * (e472 * e588)) / e378 +
                                                                                                                                                                                                                                                                                                                                                            0.5 * (e353 * e27 - (0.5 * ((e293 - 2 * (e55 / e7)) *
                                                                                                                                                                                                                                                                                                                                                                                        e63 / e7) + 0.5 * (τ * (e79 / e7 + e78) * e23 / e8)) *
                                                                                                                                                                                                                                                                                                                                                                                e28) + τ * ((0.5 * (0.5 * (e434 - e79 * e28 / e7) +
                                                                                                                                                                                                                                                                                                                                                                                                    0.5 * e434) + τ * (0.5 * (e157 * e27 - (e343 +
                                                                                                                                                                                                                                                                                                                                                                                                                                            0.5 * (e157 * e13 / e7)) * e13) - 0.5 * (e379 * e13))) / e8 +
                                                                                                                                                                                                                                                                                                                                                                                            0.5 * (e312 * e27 + 0.5 * (e301 * e31 / e8) - (((4 *
                                                                                                                                                                                                                                                                                                                                                                                                                                             (e262 * e27) - (0.375 * e98 + 1 / e7) * e28 / e7) *
                                                                                                                                                                                                                                                                                                                                                                                                                                            e18 + 0.5 * (e79 / e8)) * e13 / e7 + 2 * (e32 *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      e296)) * e13)) * e23) / e378 + (((((e482 - (0.5 *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  (e569 / e7) + e437)) * e181 - e51 * e588) / e8 - (e512 *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    e124 + 2 * (e469 * e51 * e33 * e23 * e5 / e309))) / e166 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       ((τ * ((e400 + e597 - e54 * (e314 + e255) * e13) *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              e27 - (e464 * e103 + (e119 * e345 - e54 * (0.5 *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         ((e290 + e293) / e7) + 0.5 * (τ * e116 * e23 / e8)) *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  e13) * e13 + e372 * e27 + (0.5 * e367 + 0.5 *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          (e236 + ((e251 + τ * (e227 + e529)) * e13 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   e241) / e8) + 0.5 * (τ * e119 * e13 / e8)) *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            e23 + 0.5 * (e417 / e7)) * e13) / e8 - ((((e177 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       e242 - e76) / e7 + 0.75 * e72) * e27 - ((e251 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                0.25 * (e265 / e8)) * e13 + 0.5 * e274 + 0.5 *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         ((e260 + e587) * e13 + e241) + e589) * e13 / e7) / e8 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    (e525 + (0.5 * (e71 * e31 / e8) - 2 * (e232 *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           e21)) / e20 + e530 + e531 - e537) / e20) *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   e55) * e18 - (((e328 * (e318 - e354) + (e54 *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           ((2 * (((e85 - e176) * e51 * e21 / e20 + (e150 *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     e21 + 0.5 * (e51 * e31 / e8)) * e33) / e20) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             2 * ((e573 * e23 / e20 - 0.5 * (τ * e562 / e8)) *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  e51 + e150 * (2 * e194 - e207) - e512 *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   e33) - 0.5 * (τ * e218 * e13 / e8)) * e23 *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            e74 + 0.5 * (((2 * (e92 / e20) - 0.5 * (e577 / e7)) *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          e31 - 0.5 * (e577 * e63 / e7)) / e8)) + 2 *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  (e119 * e226) - e464 * e94) * e55 + e140 *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      (0.5 * (τ * e149 * e23 / e8) - 2 * (((e226 -
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            e94^2 / e146) * e55 + 0.5 * (τ * e94 * e13 *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         e23 / e8)) / e146))) * e18 + e366 * e226 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  e94 * ((τ * ((e247 + (e529 - e487) * e23) *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               e13 + e477) / e8 - e386) * e18 - (e385 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 e329 * e32 + e594)) + 2 * (e51 * e484 *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            e33 * e23 * e74) - e364 * e32) / e146 + ((0.5 *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      ((e119 - 0.5 * (e572 / e7)) * e31 - 0.5 * (e572 *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 e63 / e7)) + 0.5 * (e119 * e31)) / e8 - e464 *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         e32) * e55 + e32 * (τ * ((e247 + e371 - e54 *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ((0.5 * e103 + 0.5 * e355) / e7 + 1.25 * e255) *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 e13) * e13 + (e119 * e13 + e54 * e27) * e355 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  e396) / e8 - 2 * e386) + (0.5 * (e422 * e31) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            0.5 * (τ * e272 * e13 * e23)) / e8)) * e5 / e59 -
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       (e507 * e21 / e20 + 0.5 * (τ * (e360 * e27 + (e560 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     ((τ * ((e247 + (e586 - e487) * e23) * e13 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            e477) / e8 - (e504 + e305) * e55) * e18 -
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      (e270 * e32 + e304 + e579)) * e5 / e59 -
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     0.5 * (e360 * e13 / e7)) * e13) / e8))) / e20 - ((((((e466 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           2 * (e374 - e87)) * e23 + e590) * e51 + e150 *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   e336) * e124 + e465 * e181 / e8) * e33 + (e194 -
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             2 * ((e166 * e21 + 2 * (e90 * e5)) * e166 * e33 *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  e20 / e424)) * e51 * e336 * e124) * e5 / e424 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      0.5 * (τ * e553 * e13 / e8))) * e23 * e5) * m.v0) *
              m.σ) * m.σ + (e547 + (e463 + ((e539 -
                                             2 * e370) / e33 - (((((e541 - (e67 + e592) * e13) / e8 +
                                                                   e601) * e55 + e604) * e18 - e599) * e5 / e59 + e256) *
                                                               e23 * e5 / e20) * m.v0 - ((e536 + 4 * e246) / e20 +
                                                                                         e593) * m.κ * m.θ) * e331 + (((e21 * e596 / e20 +
                                                                                                                        2 * (e549 + e555 - e321) + 2 * (e551 + e555 - e321)) / e20 +
                                                                                                                       e300) * m.κ * m.θ - ((e560 + ((τ * (0.5 *
                                                                                                                                                           (e407 * e13 * e23) + 2 * ((e247 - 0.5 * (e253 *
                                                                                                                                                                                                    e23 / e8)) * e13 + e477)) / e8 - ((((0.5 * e382 + e346 +
                                                                                                                                                                                                                                         e191) / e7 + e78) * e13 + e241 + e533) / e8 + 2 *
                                                                                                                                                                                                                                                                                       e504) * e55) * e18 - (e442 * e32 + 2 * (e497 *
                                                                                                                                                                                                                                                                                                                               e51 * e33 * e23 * e74 / e146) + e594)) * e5 / e59 -
                                                                                                                                             e546 * e43) * e23 * e5 / e20 + ((e51 * (e539 - e566) +
                                                                                                                                                                              e306) * e23 / e69 + 2 * e574) / e33) * m.v0) *
                                                                                                                     m.σ + e565 + ((e351 + 2 * (e543 / e20)) * m.κ *
                                                                                                                                   m.θ - (e553 * e23 * e5 + e554) * m.v0) *
                                                                                                                                  m.σ) * e331 + e503 * e474 + (e463 + ((e538 -
                                                                                                                                                                        e370) / e33 - e546 * e23 * e5 / e20) * m.v0 -
                                                                                                                                                               (e596 / e20 + e72) * m.κ * m.θ) * e524) *
     e490 + (0 - 1im) * ((e547 + e331 * (((e539 - e370) / e33 -
                                          ((((((e495 - e30) * e36 - (e66 + e592) * e13) / e8 +
                                              e601) * e55 + e604) * e18 - e599) * e5 / e59 +
                                           2 * e256) * e23 * e5 / e20) * m.v0 + 2 * e474 +
                                         4 * e463 - ((e535 + 4 * e421) / e20 + e593) * m.κ *
                                                    m.θ) + (((2 * e543 + 2 * (e500 + e551 - e321)) / e20 +
                                                             e300) * m.κ * m.θ + 2 * ((e351 + 2 * ((e549 +
                                                                                                    (e344 + 2 * e459) / e20 - e321) / e20)) * m.κ *
                                                                                      m.θ - (((e51 * (e538 - e566) + e306) * e23 / e69 +
                                                                                              e574) / e33 + ((τ * ((e247 + (0.5 * e373 - e487) *
                                                                                                                           e23) * e13 + e477) / e8 - (e504 + e427) * e55) *
                                                                                                             e18 - (e415 * (e94 / e146 + e40) + (e394 + e403) *
                                                                                                                                                e32 + e304)) * e23 * e74 / e146) * m.v0) - ((((e51 *
                                                                                                                                                                                               (e538 - e534) + 2 * e306) / e166 + ((τ * (e179 *
                                                                                                                                                                                                                                         e54 * e103 + ((e586 + e529 - e487) * e23 + e597 -
                                                                                                                                                                                                                                                       e395) * e13) / e8 - ((e461 + e21 * e598 / e20 + e139 +
                                                                                                                                                                                                                                                                             e488 + e489) / e20 + 2 * e305) * e55) * e18 - (e385 +
                                                                                                                                                                                                                                                                                                                            ((e598 / e20 + 2 * e145) * e55 + e253 * (2.5 * e23 -
                                                                                                                                                                                                                                                                                                                                                                     1) / e8) * e32 + e579 + 2 * e304)) * e5 / e59 -
                                                                                                                                                                                              e360 * (e42 + e39)) / e20 - e516) * e23 * e5 + 2 *
                                                                                                                                                                                                                                             e554) * m.v0) * m.σ + 3 * (e331 * e474 +
                                                                                                                                                                                                                                                                        e565)) * e490 + (0 - 1im) * (((e181 / e378 - e360 * e23 *
                                                                                                                                                                                                                                                                                                                     e5 / e20) * m.v0 + 5 * e474 + 6 * e463 - (e78 +
                                                                                                                                                                                                                                                                                                                                                               2 * (e421 / e20)) * m.κ * m.θ) * e490 +
                                                                                                                                                                                                                                                                                                     (0 - 1im) * ((0 - 1im) * e440 * evaluateCharFunc(cf, τ, u) + 4 * e286 * e490) * e440) * e440) * e440) *
    exp((0 - 1im) * u * e440)

end

function computeCumulants(cf::CharFunc{HestonParams{T},CR}, τ::T) where {T,CR}
    return computeCumulants(model(cf), τ)
end

function computeCumulants(p::HestonParams{T}, τ::T) where {T}
    lambda = p.κ
    ubar = p.θ
    u0 = p.v0
    eta = p.σ
    rho = p.ρ
    term = τ
    if lambda == 0 && eta == 0
        c1 = -u0 * term / 2
        c2 = u0 * term
        c4 = 0
        return c1, c2, c4
    elseif lambda == 0
        c1 = -u0 * term / 2
        c2 = u0 * term * (1 + term * term * eta * (eta * term / 12 - rho / 2))
        c4 =
            u0 *
            eta^2 *
            term^3 *
            (
                2 * rho^2 - 2 * eta * term * rho +
                eta^4 * term^4 * 17 / 1680 +
                eta^4 * term^5 * rho^2 * 11 / 20 - eta * term * rho^3 / 2 +
                eta^2 * term^2 * 3 / 10 - eta^3 * term^3 * rho * 17 / 120 + 1
            )
        return c1, c2, c4
    end
    elt = exp(-lambda * term)
    c1 = (1 - elt) * (ubar - u0) / (2 * lambda) - ubar * term / 2
    lambda2 = lambda * lambda
    lambda3 = lambda * lambda2
    eta2 = eta * eta
    elt2 = elt * elt
    # correct c2
    c2 =
        1 / (8 * lambda3) * (
            2 *
            u0 *
            (
                lambda2 * 4 * (elt * term * rho * eta + 1 - elt) +
                lambda * (4 * eta * rho * (elt - 1) - 2 * elt * term * eta2) +
                eta2 * (1 - elt2)
            ) +
            ubar * (
                8 * lambda3 * term +
                lambda2 * 8 * (-eta * rho * term * (1 + elt) + elt - 1) +
                lambda * (2 * eta2 * term * (1 + 2 * elt) + eta * rho * 16 * (1 - elt)) +
                eta2 * (-5 + 4 * elt + elt2)
            )
        )
    # paper c2
    # c2 = 1 / (8 * lambda3) * (eta*term*lambda*elt*(u0-ubar)*(8*lambda*rho-4*eta) + lambda*rho*eta*(1-elt)*(16*ubar-8*u0) + 2*ubar*lambda*term*(-4*lambda*rho*eta+eta2+4*lambda2) + eta2*((ubar-2*u0)*elt2+ubar*(6*elt-7)+2*u0) + 8*lambda2*(u0-ubar)*(1-elt))

    lambda4 = lambda2 * lambda2
    lambda5 = lambda4 * lambda
    lambda6 = lambda4 * lambda2
    lambda7 = lambda3 * lambda4
    eta3 = eta * eta2
    term2 = term * term
    term3 = term2 * term
    rho2 = rho * rho
    eta4 = eta2 * eta2
    elt3 = elt2 * elt
    elt4 = elt2 * elt2
    rho3 = rho2 * rho
    c4 =
        -1 / lambda7 *
        2 *
        ubar *
        eta2 *
        (
            (
                (term3 * rho3 * eta - 3 * term2 * rho2) * lambda6 -
                3 *
                term *
                (term2 * rho2 * eta2 - 4 * term * rho * (rho2 + 1) * eta + 8 * rho2 + 2) *
                lambda5 / 2 +
                (
                    3 * term3 * rho * eta3 / 4 - 21 * term2 * (rho2 + 3 / 14) * eta2 / 2 +
                    (18 * term * rho3 + 24 * term * rho) * eta - 18 * rho2 - 3
                ) * lambda4 -
                (
                    term3 * eta3 - 42 * term2 * rho * eta2 +
                    (240 * term * rho2 + 54 * term) * eta - 192 * rho3 - 192 * rho
                ) *
                eta *
                lambda3 / 8 -
                3 *
                eta2 *
                (term2 * eta2 - 35 * term / 2 * rho * eta + 40 * rho2 + 15 / 2) *
                lambda2 / 4 - 27 * eta3 * (term * eta - 20 * rho / 3) * lambda / 16 -
                21 * eta4 / 16
            ) * elt +
            (
                (-3 / 4 + 3 * term * rho * eta - 3 * term2 * rho2 * eta2 / 2) * lambda4 +
                3 *
                eta *
                (term2 * rho * eta2 + (-4 * term * rho2 - 3 * term / 2) * eta + 4 * rho) *
                lambda3 / 2 -
                3 *
                eta2 *
                (term2 * eta2 - 14 * term * rho * eta + 20 * rho2 + 6) *
                lambda2 / 8 + (-15 / 16 * term * eta4 + 9 * eta3 * rho / 2) * lambda -
                21 * eta4 / 32
            ) * elt2 +
            3 *
            eta2 *
            (
                (term * rho * eta - 1) * lambda2 +
                (-term / 2 * eta2 + 2 * rho * eta) * lambda - eta2 / 2
            ) *
            elt3 / 8 - 3 * eta4 * elt4 / 128 +
            (-6 * term * rho2 - 3 * term / 2) * lambda5 +
            ((6 * term * rho3 + 9 * term * rho) * eta + 18 * rho2 + 15 / 4) * lambda4 -
            9 * eta * ((rho2 + 0.25) * term * eta + 8 * rho3 / 3 + 10 * rho / 3) * lambda3 +
            15 * eta2 * (term * rho * eta + 10 * rho2 + 11 / 5) * lambda2 / 4 +
            (-33 / 2 * eta3 * rho - 15 / 32 * term * eta4) * lambda +
            279 * eta4 / 128
        )
    c4 +=
        u0 / lambda7 * (
            2 *
            eta2 *
            (
                (
                    (term3 * rho3 * eta - 3 * term2 * rho2) * lambda6 -
                    3 *
                    term *
                    (
                        term2 * rho2 * eta2 - 2 * term * rho * (rho2 + 2) * eta +
                        4 * rho2 +
                        2
                    ) *
                    lambda5 / 2 +
                    (
                        3 * term3 * rho * eta3 / 4 - 6 * (rho2 + 3 / 8) * term2 * eta2 +
                        6 * term * rho * (rho2 + 2) * eta - 6 * rho2
                    ) * lambda4 -
                    eta *
                    (
                        term3 * eta3 - 24 * term2 * rho * eta2 +
                        (72 * term * rho2 + 18 * term) * eta - 48 * rho3
                    ) *
                    lambda3 / 8 -
                    3 * eta2 * (term2 * eta2 - 7 * term * rho * eta - 3) * lambda2 / 8 -
                    3 * eta3 * (term * eta + 10 * rho) * lambda / 16 + 3 * eta4 / 8
                ) * elt +
                (
                    (6 * term * rho * eta - 3 * term2 * rho2 * eta2 - 3 / 2) * lambda4 +
                    3 *
                    (
                        term2 * rho * eta2 +
                        (-3 * term * rho2 - 3 * term / 2) * eta +
                        3 * rho
                    ) *
                    eta *
                    lambda3 -
                    3 *
                    eta2 *
                    (term2 * eta2 - 10 * term * rho * eta + 12 * rho2 + 3) *
                    lambda2 / 4 - 9 * eta3 * (term * eta - 10 * rho / 3) * lambda / 8 -
                    3 * eta4 / 8
                ) * elt2 +
                9 *
                eta2 *
                (
                    (term * rho * eta - 1) * lambda2 +
                    (-term / 2 * eta2 + 5 / 3 * rho * eta) * lambda - eta2 / 3
                ) *
                elt3 / 8 - 3 * eta4 * elt4 / 32 -
                6 *
                ((rho2 + 1 / 4) * lambda2 - 5 * lambda * rho * eta / 4 + 5 * eta2 / 16) *
                (lambda * rho * eta - eta2 / 4 - lambda2)
            )
        )
    return c1, c2, c4
end
