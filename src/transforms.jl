"""
   local_axes(m::MUonEModule)

Compute the versors of the local reference frame. 
    
"""
function local_axes(m::MUonEModule)
    # in case of wrong understanding of
    # the transform, traspose R (=R^-1)
    # and then do the same thing.
    eX = SVector{3}(m.R[1:3]) # poi si spiega
    eY = SVector{3}(m.R[4:6])
    eZ = SVector{3}(m.R[7:9])
    return eX, eY, eZ
end

"""
   local_to_global(q::StaticVector{3}, m::MUonEModule)

Passive transform of the vector `q` in the global coordinates.
 
"""
function local_to_global(q::StaticVector{3}, m::MUonEModule)
	return m.R * q + m.r0
end

"""
   global_to_local(r::StaticVector{3}, m::MUonEModule)

Passive transform of the vector `r` in the the local coordinates.
 
"""
function global_to_local(r::StaticVector{3}, m::MUonEModule)
    return inv(m.R) * (r - m.r0)
end

"""
    strip_to_local(strip_X::Real, strip_Y::Real)
    
Conversion of module local coordinates from strip number to cm.

# References
See https://gitlab.cern.ch/muesli/daq-sw/daq-decode/-/blob/api/src/Converter.cpp
"""
function strip_to_local(strip_X::Real, strip_Y::Real, m::MUonEModule)
    nstrips = 1016
    strip_pitch = 0.009
    sensor_dimension_Y = 10
    
    strip_X = strip_X/2 - 1
    
    cm_X = (nstrips/2 - strip_X) * strip_pitch - strip_pitch/2
    cm_Y = (strip_Y - 0.5) * sensor_dimension_Y
    if m.type == 'U' || m.type == 'V'
        return SVector(cm_X, cm_Y, -2.0)
    else
        return SVector(cm_X, cm_Y, -0.9)
    end
end
