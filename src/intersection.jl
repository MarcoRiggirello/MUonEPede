"""
   intersection(m::MUonEModule, t::Track)

Computes the global coordinate `z` of the intersection point between a `Track` and a `MUonEModule`.
 
"""
function intersection(m::MUonEModule, t::Track)
    _, _, ew = local_axes(m)
    Z = m.spacing/2
    return (-Z + (m.r0 - t.t0) ⋅ ew) / (t.et ⋅ew),
           (Z + (m.r0 - t.t0) ⋅ ew) / (t.et ⋅ew)
end

function interpolate(stubs::StubSet{T}, modules::MUonEStation) where T

    hit0 = local_to_global(stub_to_local(stubs[1], modules[1])[1], modules[1])
    hit1 = local_to_global(stub_to_local(stubs[2], modules[2])[1], modules[2])
    hit4 = local_to_global(stub_to_local(stubs[5], modules[5])[1], modules[5])
    hit5 = local_to_global(stub_to_local(stubs[6], modules[6])[1], modules[6])

    Δx_40 = hit4.x - hit0.x
    Δy_51 = hit5.y - hit1.y
    Δz_40 = hit4.z - hit0.z
    Δz_51 = hit5.z - hit1.z
    
    x0 = (hit4.z * hit0.x - hit0.z * hit4.x)/Δz_40 
    y0 = (hit5.z * hit1.y - hit1.z * hit5.y)/Δz_51 
     
    mx = Δx_40 / Δz_40
    my = Δy_51 / Δz_51
    
    return Track{T}(x0, y0, mx, my) 
end
