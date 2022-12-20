"""
   local_axes(m::MUonEModule)

Compute the versors of the local reference frame. 
    
"""
function local_axes(m::MUonEModule)
    eX = SVector{3}(m.R[1:3])
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
    strip_to_local(strip_X::Real, strip_Y::Real, m::MUonEModule)
    
Conversion of module local coordinates from strip number to cm.

# References
See https://gitlab.cern.ch/muesli/daq-sw/daq-decode/-/blob/api/src/Converter.cpp

# Notes
To ensure that the local frame is right handed, the cm_X output has an overall minus sign
wrt the original converter.
"""
function strip_to_local(strip_X::T, strip_Y::T, m::MUonEModule) where {T<:Real}
    nstrips = 1016
    strip_pitch = 0.009
    sensor_dimension_Y = 10

    cm_X = (strip_X - 3 - nstrips) * strip_pitch / 2 # just a refactor of original code
    cm_Y = (strip_Y - 0.5) * sensor_dimension_Y
    return SVector{3, T}(cm_X, cm_Y, -1*m.spacing/2)
end

strip_to_local(strip_X, strip_Y, m) = strip_to_local(promote(strip_X, strip_Y)..., m)
