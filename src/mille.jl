# now with the mille data creation
function rmeas(l::SVector{3, T}, z, m, t) where T
    elhit = global_to_local(t(z), m)

    rmeasX = l[1] - elhit[1]
    rmeasY = l[2] - elhit[2] 
    return rmeasX, rmeasY 
end

function sigma(s::Stub{T}, w::Real=1) where T
    pitch = 0.009

    st = @SVector [.0, .5]
    ct = @SVector [.0, .5, .25, .75]

    S = @SVector [√12, √24]
    C = @SMatrix [√12 √24 √16 √16;
                  √24 √12 √16 √16]

    sb = st .== s.localX % 1
    cb = ct .== s.bend % 1

    sf = S ⋅ sb
    cf = sb ⋅ (C * cb)
    return T(pitch * w / sf), T(1.5), T(pitch * w / cf) 
end

function derlc(z::Real, m::MUonEModule)
    Я = inv(m.R)
    
    ∂ρ_∂t0x = Я * @SVector [1, 0, 0]
    ∂ρ_∂t0y = Я * @SVector [0, 1, 0]
    ∂ρ_∂mx = z * Я * @SVector [1, 0, 0]
    ∂ρ_∂my = z * Я * @SVector [0, 1, 0]

    derlcX = [∂ρ_∂t0x[1], ∂ρ_∂t0y[1], ∂ρ_∂mx[1], ∂ρ_∂my[1]]
    derlcY = [∂ρ_∂t0x[2], ∂ρ_∂t0y[2], ∂ρ_∂mx[2], ∂ρ_∂my[2]]
    return derlcX, derlcY 
end
    
function dergl(z, m, t)
    hit = t(z)
    Sx = @SMatrix [
        0  0  0;
        0  0 -1;
        0  1  0
    ]
    Sy = @SMatrix [
        0  0  1;
        0  0  0;
       -1  0  0
    ]
    Sz = @SMatrix [
        0 -1  0;
        1  0  0;
        0  0  0
    ]
    Я = inv(m.R)

    θx, θy, θz = Rotations.params(m.R)

    Яx = RotX(-θx)
    Яy = RotY(-θy)
    Яz = RotZ(-θz)

    ∂ρ_∂x0 = -Я * @SVector [1, 0, 0]
    ∂ρ_∂y0 = -Я * @SVector [0, 1, 0]
    ∂ρ_∂z0 = -Я * @SVector [0, 0, 1]

    ∂ρ_∂θx = -Яz * Яy * Sx * Яx * (hit - m.r0)
    ∂ρ_∂θy = -Яz * Sy * Яy * Яx * (hit - m.r0)
    ∂ρ_∂θz = -Sz * Яz * Яy * Яx * (hit - m.r0)

    derglX = [∂ρ_∂x0[1], ∂ρ_∂y0[1], ∂ρ_∂z0[1], ∂ρ_∂θx[1], ∂ρ_∂θy[1], ∂ρ_∂θz[1]]
    derglY = [∂ρ_∂x0[2], ∂ρ_∂y0[2], ∂ρ_∂z0[2], ∂ρ_∂θx[2], ∂ρ_∂θy[2], ∂ρ_∂θz[2]]
    return derglX, derglY
end

function label(m::MUonEModule)
    return [i + 10*m.id for i in 101:106]
end

function mille!(glder::AbstractVector, inder::AbstractVector, s, m, t, w; cic=false)
    z_s, z_c = intersection(m, t)
    l_s, l_c = stub_to_local(s, m)
    l = cic ? 36 : 24
    o = l * m.id + 1

    sigmaX, sigmaY, sigmaB = sigma(s, w)
    
    rmeasX, rmeasY = rmeas(l_s, z_s, m, t)
    derlcX, derlcY = derlc(z_s, m)
    derglX, derglY = dergl(z_s, m, t)
    
    rmeasB, _ = rmeas(l_c, z_c, m, t)
    derlcB, _ = derlc(z_c, m)
    derglB, _ = dergl(z_c, m, t)

    # Seed layer Local X residual and derivatives
    glder[o+1] = rmeasX
    inder[o+1] = zero(Int32)

    glder[o+2:o+5] = derlcX
    inder[o+2:o+5] .= 1,2,3,4

    glder[o+6] = sigmaX
    inder[o+6] = zero(Int32)

    glder[o+7:o+12] = derglX
    inder[o+7:o+12] = label(m)

    # Correlation layer LocalX residuals and derivatives
    glder[o+13] = rmeasB
    inder[o+13] = zero(Int32)

    glder[o+14:o+17] = derlcB
    inder[o+14:o+17] .= 1,2,3,4

    glder[o+18] = sigmaB
    inder[o+18] = zero(Int32)

    glder[o+19:o+24] = derglB
    inder[o+19:o+24] = label(m)

    if cic
        # Seed layer Local Y residualis and derivatives
        glder[o+25] = rmeasY
        inder[o+25] = zero(Int32)

        glder[o+26:o+29] = derlcY
        inder[o+26:o+29] .= 1,2,3,4

        glder[o+30] = sigmaY
        inder[o+30] = zero(Int32)

        glder[o+31:o+36] = derglY
        inder[o+31:o+36] = label(m)
    end
end

