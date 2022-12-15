# now with the mille data creation
function rmeas(strip_X, strip_Y, m, t)
    mlhit = strip_to_local(strip_X, strip_Y, m)
    s = intersection(m, t)
    elhit = global_to_local(t(s), m)
    return mlhit[1] - elhit[1]
end

function derlc(m::MUonEModule, t::Track)
    s = intersection(m, t)
    sθ, cθ = sincos(t.θ) 
    sϕ, cϕ = sincos(t.ϕ) 
    T = inv(m.R)
    
    dqX_dt0x = T[1,1]
    dqX_dt0y = T[1,2]
    dq_dθ = s * T * [cθ * cϕ, cθ * sϕ, -sϕ]
    dq_dϕ = s * T * [-sθ * sϕ, sθ * cϕ, 0]
    return dqX_dt0x, dqX_dt0y, dq_dθ[1], dq_dϕ[1]
end
    
function dergl(m, t)
    s = intersection(m, t)
    hit = t(s)
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
    
    dqX_dx0 = m.R[1,1]
    dqX_dy0 = m.R[2,1]
    dqX_dz0 = m.R[3,1]
    dq_dθx = RotZ(-θz) * RotY(-θy) * Sx * RotX(-θx) * (hit - m.r0)
    dq_dθy = RotZ(-θz) * Sy * RotY(-θy) * RotX(-θx) * (hit - m.r0)
    dq_dθz = Sz * RotZ(-θz) * RotY(-θy) * RotX(-θx) * (hit - m.r0)
    return dqX_dx0, dqX_dy0, dqX_dz0, dq_dθx[1], dq_dθy[1], dq_dθz[1]
end

function label(m::MUonEModule)
    return (101, 102, 103, 104, 105, 106) .+ 10*m.id
end

function mille!(glder::AbstractVector, inder::AbstractVector, strip_X, strip_Y, m, t)
    o = 12*m.id + 1

    glder[o+1] = rmeas(strip_X, strip_Y, m, t)
    inder[o+1] = zero(Int32)

    glder[o+2:o+5] .= derlc(m, t)
    inder[o+2:o+5] .= 1,2,3,4

    glder[o+6] = 0.0012 # sigma -- TO BE CHECKED
    inder[o+6] = zero(Int32)

    glder[o+7:o+12] .= dergl(m, t)
    inder[o+7:o+12] .= label(m)
end

