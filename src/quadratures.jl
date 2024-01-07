## Hermite quintic quadrature for Flinn


const p2 = (-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0)
const h0 = (
    0.0762880721720079,
    -0.0329193509714785,
    -0.3438600114643695,
    0.5060431337241802,
    1.5888963130793203,
    0.5060431337241802,
    -0.3438600114643695,
    -0.0329193509714785,
    0.0762880721720079,
)
const h1 = (
    1.60554439048059e-03,
    -3.96727024402095e-02,
    -1.57110962304862e-01,
    -2.32981551813796e-01,
    0.0,
    2.32981551813796e-01,
    1.57110962304862e-01,
    3.96727024402095e-02,
    -1.60554439048059e-03,
)

function integrateQuinticHermite(
    f,
    a::T,
    b::T,
    tol::T,
    maxRecursionDepth::Int;
    integralEstimate::T=Base.zero(T)
)::T where {T}
    if tol < eps(T)
        tol = eps(T)
    end

    if b <= a
        return Base.zero(T)
    end
    c = (a + b) / 2
    h = (b - a) / 2

    is = Base.zero(T)
    x = tuplejoin(a, c .+ h .* p2[2:end-1], b)
    y = Vector{T}(undef, length(p2))
    yp = Vector{T}(undef, length(p2))
    @inbounds for i = 1:length(p2)
        y[i], yp[i] = f(x[i])
        is += h * (h0[i] * y[i] + h * h1[i] * yp[i])
    end
    if integralEstimate == Base.zero(T)
        if is == Base.zero(T)
            is = one(T)
        end
        is = is * tol / eps(T)
    else
        is = integralEstimate * tol / eps(T)
    end

    @inbounds result =
        adaptiveQuinticHermiteAux(
            f,
            x[1],
            x[3],
            is,
            y[1],
            y[3],
            y[2],
            yp[1],
            yp[3],
            yp[2],
            maxRecursionDepth,
        ) +
        adaptiveQuinticHermiteAux(
            f,
            x[3],
            x[5],
            is,
            y[3],
            y[5],
            y[4],
            yp[3],
            yp[5],
            yp[4],
            maxRecursionDepth,
        ) +
        adaptiveQuinticHermiteAux(
            f,
            x[5],
            x[7],
            is,
            y[5],
            y[7],
            y[6],
            yp[5],
            yp[7],
            yp[6],
            maxRecursionDepth,
        ) +
        adaptiveQuinticHermiteAux(
            f,
            x[7],
            x[9],
            is,
            y[7],
            y[9],
            y[8],
            yp[7],
            yp[9],
            yp[8],
            maxRecursionDepth,
        )
    return result
end

function adaptiveQuinticHermiteAux(
    f,
    a::T,
    b::T,
    is::T,
    fa::T,
    fb::T,
    fc::T,
    fpa::T,
    fpb::T,
    fpc::T,
    bottom::Int,
)::T where {T}
    c = (a + b) / 2
    h = b - a
    d = (a + c) / 2
    e = (c + b) / 2
    fd, fpd = f(d)
    fe, fpe = f(e)
    i1 = (h / 30) * (7 * fa + 16 * fc + 7 * fb + h / 2 * (fpa - fpb))
    Sleft = (h / 60) * (7 * fa + 16 * fd + 7 * fc + h / 4 * (fpa - fpc))
    Sright = (h / 60) * (7 * fc + 16 * fe + 7 * fb + h / 4 * (fpc - fpb))
    i2 = Sleft + Sright
    if ((is + (i1 - i2) / 63) == is) || (d == e)
        return i2
    end
    if bottom <= 0
        # println("Number of maximum recursions reached")
        return i2
    end
    result =
        adaptiveQuinticHermiteAux(f, a, c, is, fa, fc, fd, fpa, fpc, fpd, bottom - 1) +
        adaptiveQuinticHermiteAux(f, c, b, is, fc, fb, fe, fpc, fpb, fpe, bottom - 1)
    return result
end

function quinticHermiteAux(
    a::T,
    b::T,
    fa::T,
    fb::T,
    fc::T,
    fpa::T,
    fpb::T,
    fpc::T,
)::T where {T}
    c = (a + b) / 2
    h = b - a
    i1 = (h / 30) * (7 * fa + 16 * fc + 7 * fb + h / 2 * (fpa - fpb))
    return i1
end


## Simpson quadratures for Filon




const p2_modsim = (-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0) ./ 4
const w2_modsim = (989.0, 5888.0, -928.0, 10496.0, -4540.0, 10496.0, -928.0, 5888.0, 989.0) ./ 14175

const p1_modsim = (-1.0, -0.5, 0.0, 0.5, 1.0)
const w1_modsim = (7.0, 32.0, 12.0, 32.0, 7.0) ./ 45
const nw1_modsim = [0.1271031188518578 -0.5084124754074312 0.7626187131111468 -0.5084124754074312 0.1271031188518578
    0.3362832433427012 -0.6725664866854024 0.0 0.6725664866854024 -0.3362832433427012
    0.5684224278099782 -0.2842112139049893 -0.5684224278099774 -0.2842112139049893 0.5684224278099782
    0.6725664866854023 0.3362832433427017 0.0 -0.3362832433427017 -0.6725664866854023]
function modsim(f, a::T, b::T, tol::T; maxRecursionDepth::Int=20)::T where {T}
    if tol < eps(T)
        tol = sqrt(eps(T))
    end
    if b <= a
        return Base.zero(T)
    end
    h = (b - a) / 2
    c = (b + a) / 2
    x = @. c + h * p2_modsim
    y = f.(x)
    is = Base.zero(T)
    isabs = Base.zero(T)
    @inbounds for i = 1:length(w2_modsim)
        is += h * w2_modsim[i] * y[i]
        isabs += h * abs(w2_modsim[i]) * abs(y[i])
    end
    noise = 50 * eps(T) * isabs
    if abs(is) <= noise
        is = b - a
    end
    tol = max(tol, noise / abs(is))

    Q = simpr(f, x[1], x[3], y[1], y[2], y[3], is, tol, maxRecursionDepth) + simpr(f, x[3], x[5], y[3], y[4], y[5], is, tol, maxRecursionDepth) + simpr(f, x[5], x[7], y[5], y[6], y[7], is, tol, maxRecursionDepth) + simpr(f, x[7], x[9], y[7], y[8], y[9], is, tol, maxRecursionDepth)
    return Q
end


function simpr(
    f,
    a::T,
    b::T,
    fa::T,
    fc::T,
    fb::T,
    is::T,
    tol::T,
    bottom::Int)::T where {T}
    h = (b - a) / 2
    c = (a + b) / 2
    x = @. c + p1_modsim * h
    fz1 = f(x[2])
    fz2 = f(x[4])
    y = (fa, fz1, fc, fz2, fb)
    fcl = y[2]
    fcr = y[4]
    i = Base.zero(T)
    iabs = Base.zero(T)
    @inbounds for k = 1:length(w2_modsim)
        i += h * w1_modsim[k] * y[k]
        iabs += h * abs(w1_modsim[k]) * abs(y[k])
    end
    noise = 50 * iabs * eps(T)
    E = h * abs.(nw1_modsim * collect(y))

    Emin = minimum(E[2:4])
    rmax = 2 * one(T)
    if iszero(Emin)
        rmax = 2.0
    else
        r = E[1:3] ./ E[2:4]
        if sum(isinf.(r)) > 0
            rmax = 2 * one(T)
        else
            rmax = maximum(r)
        end
    end
    C = 32
    err = zero(T)
    if rmax > one(T)
        err = C * maximum(E)
    elseif one(T) / 2 < rmax
        err = C * rmax * E[2]
    else
        err = C * (2 * rmax)^2 * rmax * E[2]
    end

    if E[1] < noise && E[2] < noise
        err = zero(T)
    end
    if err <= max(tol * abs(is), noise) || x[2] <= a || b <= x[4] || bottom == 0 || isnan(is) || isinf(is)
        # if ((x[2] <= a) || (b <= x[4])) && (termination2 == 0)
        #     #warning(['Interval too small: required tolerance may not be met.'])
        #     #termination2 = 1
        # end
        Q = i
    else
        Q = simpr(f, a, c, fa, fcl, fc, is, tol, bottom - 1) + simpr(f, c, b, fc, fcr, fb, is, tol, bottom - 1)
    end
    return Q
end



function integrateSimpsonGG(f, a::T, b::T, tol::T; maxRecursionDepth::Int=20, integralEstimate=zero(T))::T where {T}
    if b <= a || isnan(tol) || isinf(tol)
        return Base.zero(T)
    end

    if tol == zero(T)
        tol= eps(T)
    end
    c = (a + b) / 2

    fa = f(a)
    fm = f(c)
    fb = f(b)
    is = if (integralEstimate == zero(T))
        yy0 = f(a + 0.9501 * (b - a))
        yy1 = f(a + 0.2311 * (b - a))
        yy2 = f(a + 0.6068 * (b - a))
        yy3 = f(a + 0.4860 * (b - a))
        yy4 = f(a + 0.8913 * (b - a))

        isLocal = (b - a) / 8 * (fa + fm + fb + yy0 + yy1 + yy2 + yy3 + yy4)
        if isLocal == 0
            isLocal = b - a
        end
        isLocal
    else
        integralEstimate
    end
    is = is * tol / eps(T)

    result = adaptiveSimpsonAux(f, a, b, fa, fm, fb, is, maxRecursionDepth)
    return result
end

function adaptiveSimpsonAux(
    f,
    a::T,
    b::T,
    fa::T,
    fm::T,
    fb::T,
    is::T,
    bottom::Int,
)::T where {T}
    m = (b + a) / 2
    h = (b - a) / 4
    fml = f(a + h)
    fmr = f(b - h)
    i1 = h / 1.5 * (fa + 4 * fm + fb)
    i2 = h / 3 * (fa + 4 * (fml + fmr) + 2 * fm + fb)
    i1 = (16 * i2 - i1) / 15
    if ((is + (i1 - i2)) == is) || (m <= a) || (b <= m) || isnan(i1) || isinf(i1)
        if ((m <= a) || (b <= m))
            #Interval contains no more machine numbers, required tolerance may not be met
        end
        return i2
    end

    if bottom <= 0
        # println("Number of maximum recursions reached")
        return i2
    end

    return adaptiveSimpsonAux(f, a, m, fa, fml, fm, is, bottom - 1) +
           adaptiveSimpsonAux(f, m, b, fm, fmr, fb, is, bottom - 1)
end


function simpsonAux(
    a::T,
    b::T,
    fa::T,
    fb::T,
    fm::T
)::T where {T}
    (b - a) / 6 * (fa + 4 * fm + fb)
end



const nulw_coteda = [0.1101854769234527e-01 -0.8814838153876216e-01 0.3085193353856676 -0.6170386707713351 0.7712983384641688
    0.4267465171184540e-01 -0.2560479102710724 0.5974451239658356 -0.5974451239658356 0.0
    0.1123675793894425 -0.4775622124051307 0.6180216866419335 0.2809189484736109e-01 -0.5618378969472148
    0.2311270062743306 -0.6355992672544093 0.2311270062743311 0.5200357641172426 0.0
    0.3911196434508170 -0.5866794651762256 -0.3073082912827879 0.2514340565041003 0.5028681130081833
    0.5561911416025476 -0.2780955708012750 -0.5164632029166498 -0.3575514481730693 0.0
    0.6647755647017197 0.1661938911754310 -0.1899358756290700 -0.4036137357117497 -0.4748396890726778
    0.6455025992424890 0.4841269494318643 0.3227512996212468 0.1613756498106190 0.0]

const nw2_coteda = [nulw_coteda[1:8, 1:5] nulw_coteda[1:8, 4:-1:1]]
nw2_coteda[2:2:8, 6:9] = -nw2_coteda[2:2:8, 6:9]


function coteda(f, a::T, b::T, tol::T; maxRecursionDepth::Int=20)::T where {T}
    if tol < eps(T)
        tol = sqrt(eps(T))
    end
    if b <= a
        return Base.zero(T)
    end
    h = (b - a) / 2
    c = (b + a) / 2
    x = @. c + h * p2_modsim
    y = @. f(x)
    is = Base.zero(T)
    isabs = Base.zero(T)
    @inbounds for i = 1:length(w2_modsim)
        is += h * w2_modsim[i] * y[i]
        isabs += h * abs(w2_modsim[i]) * abs(y[i])
    end
    noise = 50 * eps(T) * isabs
    if abs(is) <= noise
        is = b - a
    end
    tol = max(tol, noise / abs(is))

    Q = coter(f, a, c, collect(y[1:5]), is, tol, maxRecursionDepth) + coter(f, c, b, collect(y[5:9]), is, tol, maxRecursionDepth)
    return Q
end


function coter(
    f,
    a::T,
    b::T,
    y::AbstractArray{T},
    is::T,
    tol::T,
    bottom::Int)::T where {T}
    h = (b - a) / 2
    c = (a + b) / 2
    x = @. c + p1_modsim * h

    i = Base.zero(T)
    iabs = Base.zero(T)
    @inbounds for k = 1:length(w1_modsim)
        i += h * w1_modsim[k] * y[k]
        iabs += h * abs(w1_modsim[k]) * abs(y[k])
    end
    noise = 50 * iabs * eps(T)
    E = h * abs.(nw1_modsim * collect(y))

    Emin = minimum(E[2:4])
    rmax = 2 * one(T)
    if iszero(Emin)
        rmax = 2.0
    else
        r = E[1:3] ./ E[2:4]
        if sum(isinf.(r)) > 0
            rmax = 2 * one(T)
        else
            rmax = maximum(r)
        end
    end
    C = 32
    err = zero(T)
    if rmax > one(T)
        err = C * maximum(E)
    elseif one(T) / 2 < rmax
        err = C * rmax * E[2]
    else
        err = C * (2 * rmax)^2 * rmax * E[2]
    end

    if E[1] < noise && E[2] < noise
        err = zero(T)
    end
    if err <= max(tol * abs(is), noise) || x[2] <= a || b <= x[4] || bottom == 0 || isnan(is) || isinf(is)
        # if ((x[2] <= a) || (b <= x[4])) && (termination2 == 0)
        #     #warning(['Interval too small: required tolerance may not be met.'])
        #     #termination2 = 1
        # end
        Q = i
    else

        z = @. c + (-3.0, -1.0, 1.0, 3.0) * h / 4
        fz = @. f.(z)
        y = (y[1], fz[1], y[2], fz[2], y[3], fz[3], y[4], fz[4], y[5])
        i = Base.zero(T)
        iabs = Base.zero(T)
        @inbounds for k = 1:length(w2_modsim)
            i += h * w2_modsim[k] * y[k]
            iabs += h * abs(w2_modsim[k]) * abs(y[k])
        end
        noise = 50 * iabs * eps(T)
        e = h * abs.(nw2_coteda * collect(y))
        e2 = e .^ 2
        E2 = @. (e2[1:2:7] + e2[2:2:8])
        E = sqrt.(E2)
        Emin = minimum(E[2:4])
        rmax = 2 * one(T)
        if iszero(Emin)
            rmax = 2.0
        else
            r = E[1:3] ./ E[2:4]
            if sum(isinf.(r)) > 0
                rmax = 2 * one(T)
            else
                rmax = maximum(r)
            end
        end
        D = 128
        err = zero(T)
        if rmax > one(T)
            err = D * maximum(E)
        elseif one(T) / 4 < rmax
            err = D * rmax * E[1]
        else
            err = D * (4 * rmax) * rmax * E[1]
        end

        if E[1] < noise && E[2] < noise
            err = zero(T)
        end
        if err <= max(tol * abs(is), noise) || z[1] <= a || b <= z[4] || bottom == 0
            # if ((x[2] <= a) || (b <= x[4])) && (termination2 == 0)
            #     #warning(['Interval too small: required tolerance may not be met.'])
            #     #termination2 = 1
            # end
            Q = i
        else
            Q = coter(f, a, c, collect(y[1:5]), is, tol, bottom - 1) + coter(f, c, b, collect(y[5:9]), is, tol, bottom - 1)
        end
    end
    return Q
end
