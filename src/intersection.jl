"""
   intersection(m::MUonEModule, t::Track)

Computes the arc lenght `s` of the intersection point between a `Track` and a `MUonEModule`.
 
"""
function intersection(m::MUonEModule, t::Track)
    _, _, eZ = local_axes(m)
    if m.type == 'X' || m.type == 'Y'
        Z = -0.9
    else
        Z = -2.0
    end
    return (Z + (m.r0 - t.t0) ⋅ eZ) / (t.et ⋅eZ)
end

function interpolate(hit0::T, hit1::T, hit4::T, hit5::T) where {T<:StaticVector{3}}
    Δx = hit4.x - hit0.x
    Δy = hit5.y - hit1.y
    Δz_40 = hit4.z - hit0.z
    Δz_51 = hit5.z - hit1.z
    
    x0 = (hit4.z * hit0.x - hit0.z * hit4.x)/Δz_40 
    y0 = (hit5.z * hit1.y - hit1.z * hit5.y)/Δz_51 
     
    xx = Δx * Δz_51
    yy = Δy * Δz_40
    zz = Δz_40 * Δz_51
    
    θ = acos(zz/hypot(xx, yy, zz))
    ϕ = atan(yy, xx)
    
    return Track(x0, y0, θ, ϕ) 
end
