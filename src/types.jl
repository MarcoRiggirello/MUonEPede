"""
    struct MUonEModule{T<:Real}
    
A 2S module used in the beam test station.
"""
struct MUonEModule{T<:Real}
    """
       id::Integer
    
    The module linkID. 
    """
    id::Integer
    """
        name::String
        
    The module name.
    """
    name::String
    """
        type::Char
        
    The module type. It must be 'X', 'Y', 'U' or 'V'.
    """
    type::Char
    """
        r0::SVector{3,T}
        
    The offset vector of the module center w.r.t. the global axes origin, in cm.
    """
    r0::SVector{3, T}
    """
        R::Rotation{3, T}
        
    The rotation matrix to be applied for local-to-global coordinate conversion.
    We use as parameters the
    [Tait-Bryan angles](https://en.wikipedia.org/wiki/Euler_angles#Tait%E2%80%93Bryan_angles)
    applied in the (right to left) order x-y-z.
    """ 
    R::Rotation{3, T}
    """
        MUonEModule{T}(x0, y0, z0, θx, θy, θz; id, name, type) where {T}
        
    default constructor
    
    #Arguments
    - `x0::Real`: the x component of `r0`, in cm;
    - `y0::Real`: the y component of `r0`, in cm;
    - `z0::Real`: the z component of `r0`, in cm;
    - `θx::Real`: rotation around global x axis, in radiants;
    - `θy::Real`: rotation around global y axis, in radiants;
    - `θz::Real`: rotation around global z axis, in radiants;
    - `id::Integer`: the linkId;
    - `name::String`: module name;
    - `type::Char`: module type. Must be 'X', 'Y', 'U', 'V',
    otherwise throws an error.
    
    #Notes
    If the module is of type 'Y' a rotation of `pi` about the `y` axis
    and a rotation of `pi/2` about the `z` axis is added, to adapt the
    constructor to the conventions used in the MUonE structure definition.
    """
    function MUonEModule{T}(x0, y0, z0, θx, θy, θz; id, name, type) where {T}
        if type == 'Y'
            θx = θx
            θy = θy + pi
            θz = θz + pi/2
        elseif type == 'X' || type == 'U' || type == 'V'
            θx = θx
            θy = θy
            θz = θz
        else
            throw(ArgumentError("module type not known."))
        end
        return new{T}(id, name, type, SVector{3, T}(x0, y0, z0), RotXYZ{T}(θx, θy, θz))
    end
end

MUonEModule(x0::T, y0::T, z0::T, θx::T, θy::T, θz::T; id, name, type) where {T<:Real} = MUonEModule{T}(x0, y0, z0, θx, θy, θz; id=id, name=name, type=type)

MUonEModule(x0::Real, y0::Real, z0::Real, θx::Real, θy::Real, θz::Real; id, name, type) = MUonEModule(promote(x0, y0, z0, θx, θy, θz)...; id=id, name=name, type=type)


"""
    struct Track{T<:Real}
    
The track model for the 160GeV muon beam.
"""
struct Track{T<:Real}
    """
        θ::T
    
    The polar angle (`z` is the polar axis), in radians.
    """
    θ::T
    """
        ϕ::T
    
    The azymuthal angle, in radians.
    """
    ϕ::T
    """
        t0::SVector{3, T}
    
    Intersection point between the track and the plane `z=0`.
    """
    t0::SVector{3, T}
    """
        et::SVector{3, T}
    
    Versor of the track direction.
    """
    et::SVector{3, T}
    """
        Track{T}(x0, y0, θ, ϕ) where {T}
    
    Default constructor.
    
    #Arguments
    - `x0::Real`: x component of `t0`, in cm;
    - `y0::Real`: y component of `t0`, in cm;
    - `θ::Real`: polar angle, in radians;
    - `ϕ::Real`: azymuthal angle, in radians.
    """
    function Track{T}(x0, y0, θ, ϕ) where {T}
        sθ, cθ = sincos(θ) 
        sϕ, cϕ = sincos(ϕ) 
        return new{T}(θ, ϕ, SVector{3, T}(x0, y0, zero(T)), SVector{3, T}(sθ*cϕ, sθ*sϕ, cθ))
    end
end

Track(x0::T, y0::T, θ::T, ϕ::T) where {T<:Real} = Track{T}(x0, y0, θ, ϕ)

Track(x0::Real, y0::Real, θ::Real, ϕ::Real) = Track(promote(x0, y0, θ, ϕ)...)

"""
    (t::Track)(s::Real)

computes the track point at a given arc length `s`.
"""
function (t::Track)(s::Real)
    return t.t0 + s * t.et
end
