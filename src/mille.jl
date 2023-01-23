# now with the mille data creation
function rmeas(s, z, m, t)
    mlhit = strip_to_local(s, m)
    elhit = global_to_local(t(z), m)
    
    rmeasX = mlhit[1] - elhit[1]
    rmeasY = mlhit[2] - elhit[2] 
    return rmeasX, rmeasY 
end

function sigma(w)
    return Float32(0.01f0*w), 1.5f0 
end

function derlc(z::Real, m::MUonEModule)
    dqX_dt0x = m.R[1,1]
    dqX_dt0y = m.R[2,1]
    dqX_dmx = z * m.R[1,1]
    dqX_dmy = z * m.R[2,1]

    dqY_dt0x = m.R[1,2]
    dqY_dt0y = m.R[2,2]
    dqY_dmx = z * m.R[1,2]
    dqY_dmy = z * m.R[2,2]

    derlcX = [dqX_dt0x, dqX_dt0y, dqX_dmx, dqX_dmy]
    derlcY = [dqY_dt0x, dqY_dt0y, dqY_dmx, dqY_dmy]
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
    θx, θy, θz = Rotations.params(m.R)

    Rx = RotX(-θx)
    Ry = RotY(-θy)
    Rz = RotZ(-θz)

    dqX_dx0 = m.R[1,1]
    dqX_dy0 = m.R[2,1]
    dqX_dz0 = m.R[3,1]

    dqY_dx0 = m.R[1,2]
    dqY_dy0 = m.R[2,2]
    dqY_dz0 = m.R[3,2]

    dq_dθx = Rz * Ry * Sx * Rx * (hit - m.r0)
    dq_dθy = Rz * Sy * Ry * Rx * (hit - m.r0)
    dq_dθz = Sz * Rz * Ry * Rx * (hit - m.r0)

    derglX = [dqX_dx0, dqX_dy0, dqX_dz0, dq_dθx[1], dq_dθy[1], dq_dθz[1]]
    derglY = [dqY_dx0, dqY_dy0, dqY_dz0, dq_dθx[2], dq_dθy[2], dq_dθz[2]]
    return derglX, derglY
end

function label(m::MUonEModule)
    return [i + 10*m.id for i in 101:106]
end

function mille!(glder::AbstractVector, inder::AbstractVector, s, m, t, w)
    z = intersection(m, t)
    o = 24*m.id + 1

    rmeasX, rmeasY = rmeas(s, z, m, t)
    sigmaX, sigmaY = sigma(w)
    derlcX, derlcY = derlc(z, m)
    derglX, derglY = dergl(z, m, t)

    # Local X residual and derivatives
    glder[o+1] = rmeasX
    inder[o+1] = zero(Int32)

    glder[o+2:o+5] = derlcX
    inder[o+2:o+5] .= 1,2,3,4

    glder[o+6] = sigmaX
    inder[o+6] = zero(Int32)

    glder[o+7:o+12] = derglX
    inder[o+7:o+12] = label(m)

    # Local Y residual and derivatives
    glder[o+13] = rmeasY
    inder[o+13] = zero(Int32)

    glder[o+14:o+17] = derlcY
    inder[o+14:o+17] .= 1,2,3,4

    glder[o+18] = sigmaY
    inder[o+18] = zero(Int32)

    glder[o+19:o+24] = derglY
    inder[o+19:o+24] = label(m)
end

