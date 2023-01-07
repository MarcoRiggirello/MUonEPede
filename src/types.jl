"""

    MUonEModule{T<:Real}(x0, y0, z0, θx, θy, θz; id, name, spacing)
    
A 2S module used in the beam test station.

`x0`, `y0` and `z0` are used to construct the (static) vector `r0`,
which defines the position of the center of the module in the global
reference frame.

The `θ`s are the [Tait-Bryan angles](https://en.wikipedia.org/wiki/Euler_angles#Tait%E2%80%93Bryan_angles)
which defines the active, extrinsic, local to global rotation
matrix  `R`, applied in the (right to left) order x-y-z.

#Arguments
- `x0::T`: the x component of `r0`, in cm;
- `y0::T`: the y component of `r0`, in cm;
- `z0::T`: the z component of `r0`, in cm;
- `θx::T`: rotation around global x axis, in radiants;
- `θy::T`: rotation around global y axis, in radiants;
- `θz::T`: rotation around global z axis, in radiants;
- `id::Integer`: the linkId;
- `name::String`: module name;
- `spacing::T`: spacing between sensors in cm.
"""
struct MUonEModule{T<:Real}
    id::Integer
    name::String
    spacing::T
    r0::SVector{3, T}
    R::Rotation{3, T}
    function MUonEModule{T}(x0, y0, z0, θx, θy, θz; id, name, spacing) where {T}
        return new{T}(id, name, spacing, SVector{3, T}(x0, y0, z0), RotXYZ{T}(θx, θy, θz))
    end
end

MUonEModule(x0::T, y0::T, z0::T, θx::T, θy::T, θz::T; id, name, spacing::T) where {T<:Real} = MUonEModule{T}(x0, y0, z0, θx, θy, θz; id=id, name=name, spacing=spacing)

function MUonEModule(x0::Real, y0::Real, z0::Real, θx::Real, θy::Real, θz::Real; id, name, spacing::Real)
    x0, y0, z0, θx, θy, θz, spacing = promote(x0, y0, z0, θx, θy, θz, spacing)
    return MUonEModule(x0, y0, z0, θx, θy, θz; id=id, name=name, spacing=spacing)
end

"""

    MUonEStation{T}(x1, y1, u, v, x2, y2)
    MUonEStation(x1, y1, u, v, x2, y2)

A StaticArrays.FieldVector of MUonEModules defining a station.
"""
struct MUonEStation{T} <: FieldVector{6, MUonEModule{T}}
    x1::MUonEModule{T}
    y1::MUonEModule{T}
    u::MUonEModule{T}
    v::MUonEModule{T}
    x2::MUonEModule{T}
    y2::MUonEModule{T}
end

"""

    Stub{T}(localX, localY, bend, link)
    Stub(localX, localY, bend, link)
    Stub(;bend, link, localX, localY)

The stub data structure.

#Arguments
to be completed
"""
struct Stub{T<:Real}
    localX::T
    localY::T
    bend::T
    link::Integer
end

Stub(localX::Real, localY::Real, bend::Real, link) = Stub(promote(localX, localY, bend)..., link)

Stub(;bend, link, localX, localY) = Stub(localX, localY, bend, link)

"""

    Track{T<:Real}(x0::T, y0::T, mx::T, my::T)
    
The track model for the 160GeV muon beam.
"""
struct Track{T<:Real}
    t0::SVector{3, T}
    et::SVector{3, T}
    function Track{T}(x0, y0, mx, my) where {T}
        return new{T}(SVector{3, T}(x0, y0, zero(T)), SVector{3, T}(mx, my, one(T)))
    end
end

Track(x0::T, y0::T, mx::T, my::T) where {T<:Real} = Track{T}(x0, y0, mx, my)

Track(x0::Real, y0::Real, mx::Real, my::Real) = Track(promote(x0, y0, mx, my)...)

"""

    (t::Track)(z::Real)

computes the track point at a given global coordinate `z`.
"""
function (t::Track)(z::Real)
    return t.t0 + z * t.et
end

