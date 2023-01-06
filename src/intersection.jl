"""
   intersection(m::MUonEModule, t::Track)

Computes the global coordinate `z` of the intersection point between a `Track` and a `MUonEModule`.
 
"""
function intersection(m::MUonEModule, t::Track)
    _, _, eZ = local_axes(m)
    Z = -1 * m.spacing/2
    return (Z + (m.r0 - t.t0) ⋅ eZ) / (t.et ⋅eZ)
end

function interpolate(stubs::AbstractVector{Stub}, modules::AbstractVector{MUonEModule})

    hit0 = local_to_global(strip_to_local(stubs[1], modules[1]), modules[1])
    hit1 = local_to_global(strip_to_local(stubs[2], modules[2]), modules[2])
    hit4 = local_to_global(strip_to_local(stubs[5], modules[5]), modules[5])
    hit5 = local_to_global(strip_to_local(stubs[6], modules[6]), modules[6])

    Δx_40 = hit4.x - hit0.x
    Δy_51 = hit5.y - hit1.y
    Δz_40 = hit4.z - hit0.z
    Δz_51 = hit5.z - hit1.z
    
    x0 = (hit4.z * hit0.x - hit0.z * hit4.x)/Δz_40 
    y0 = (hit5.z * hit1.y - hit1.z * hit5.y)/Δz_51 
     
    mx = Δx_40 / Δz_40
    my = Δy_51 / Δz_51
    
    return Track(x0, y0, mx, my) 
end
