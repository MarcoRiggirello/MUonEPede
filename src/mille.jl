# now with the mille data creation
function rmeas(l::SVector{3, T}, z, m, t) where T
    elhit = global_to_local(t(z), m)

    rmeasX = l[1] - elhit[1]
    rmeasY = l[2] - elhit[2] 
    return rmeasX, rmeasY 
end

function sigma(w)
    return Float32(0.003*w), 1.5f0, Float32(0.01*w) 
end

function derlc(z::Real, m::MUonEModule)
    Я = inv(m.R)
    
    dρ_dt0x = Я * @SVector [1, 0, 0]
    dρ_dt0y = Я * @SVector [0, 1, 0]
    dρ_dmx = z * Я * @SVector [1, 0, 0]
    dρ_dmy = z * Я * @SVector [0, 1, 0]

    derlcX = [dρ_dt0x[1], dρ_dt0y[1], dρ_dmx[1], dρ_dmy[1]]
    derlcY = [dρ_dt0x[2], dρ_dt0y[2], dρ_dmx[2], dρ_dmy[2]]
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

    Rx = RotX(θx)
    Ry = RotY(θy)
    Rz = RotZ(θz)

    dρ_dx0 = -Я * @SVector [1, 0, 0]
    dρ_dy0 = -Я * @SVector [0, 1, 0]
    dρ_dz0 = -Я * @SVector [0, 0, 1]

    dρ_dθx = -Я * (Sx * Rx * Ry * Rz) * Я  * (hit - m.r0)
    dρ_dθy = -Я * (Rx * Sy * Ry * Rz) * Я  * (hit - m.r0)
    dρ_dθz = -Я * (Rx * Ry * Sz * Rz) * Я  * (hit - m.r0)

    # dρ_dθx = Я * (Sx * m.R) * Я * (hit - m.r0)
    # dρ_dθy = Я * (Sy * m.R) * Я * (hit - m.r0)
    # dρ_dθz = Я * (Sz * m.R) * Я * (hit - m.r0)

    derglX = [dρ_dx0[1], dρ_dy0[1], dρ_dz0[1], dρ_dθx[1], dρ_dθy[1], dρ_dθz[1]]
    derglY = [dρ_dx0[2], dρ_dy0[2], dρ_dz0[2], dρ_dθx[2], dρ_dθy[2], dρ_dθz[2]]
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

    sigmaX, sigmaY, sigmaB = sigma(w)
    
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

