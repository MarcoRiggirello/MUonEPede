# a dummy Monte Carlo for pede alignment

"""
   Do the inverse of strip_to_local 
"""
function local_to_strip(q::StaticVector{3,T}) where {T<:Real}
    nstrips = 1016
    strip_pitch = 0.009

    strip_X = 2 * q.x ÷ strip_pitch + nstrips + 3
    strip_Y = q.y > 0 ? 0.75 : 0.25
    return convert.(T, [strip_X, strip_Y])
end

function mcdata!(strips_X, strips_Y, modules::StaticVector{6, MUonEModule})
    θ = 1.2e-3 * rand(Float32)
    ϕ = pi * rand(Float32)
    x0 = randn(Float32)
    y0 = randn(Float32)
    
    t = Track{Float32}(x0, y0, θ, ϕ)
    s = intersection.(modules, [t])
    
    vov = local_to_strip.(global_to_local.(t.(s), modules))
    strips_X .= vov[1][1]
    strips_Y .= vov[1][2]
end

"""
    generatebinmc(; nevents::Integer, mcfname::String, nmfname::String, ofname::String)

Generate MonteCarlo data in a Fortran binary file using 
as input the XML structure file.

#Arguments
- `nevents::Integer`: number of tracks;
- `mcfname::String`: name of the xml structure file with the montecarlo truth;
- `nmfnmae::String`: name of the xml structure file with nominal positions; 
- `ofname::String`: name of the fortran binary output file.
"""
function generatebinmc(; nevents::Integer, mcfname::String, nmfname::String, ofname::String)
    ffile = FortranFile(ofname, "w")
    mcmodules = getmodules(mcfname)
    nmmodules = getmodules(nmfname)
    
    # service arrays
    strips_X = MVector{6, Float32}(undef)
    strips_Y = MVector{6, Float32}(undef)

    # see https://gitlab.desy.de/claus.kleinwort/millepede-ii/-/blob/main/mille.f90#L27
    S = 73 # (1 rmeas + 1 sigma + 4 lder + 6 gder) * 6 measurements + 1 line of zeros
    glder = Vector{Float32}(undef, S)
    inder = Vector{Int32}(undef, S)

    glder[1] = zero(Float32)
    inder[1] = zero(Int32)
    for _ in ProgressBar(1:nevents)
        mcdata!(strips_X, strips_Y, mcmodules)
        # costruisci la traccia target
        hit0 = local_to_global(strip_to_local(strips_X[1], strips_Y[1], nmmodules[1]), nmmodules[1])
        hit1 = local_to_global(strip_to_local(strips_X[2], strips_Y[2], nmmodules[2]), nmmodules[2])
        hit4 = local_to_global(strip_to_local(strips_X[5], strips_Y[5], nmmodules[5]), nmmodules[5])
        hit5 = local_to_global(strip_to_local(strips_X[6], strips_Y[6], nmmodules[6]), nmmodules[6])
        track = interpolate(hit0, hit1, hit4, hit5)
        # mille() per ogni hit
        for (sx, sy, m) in zip(strips_X, strips_Y, nmmodules)
            mille!(glder, inder, sx, sy, m, track)
        end
        # scrivi su file
        write(ffile, Int32(S+S), glder, inder) # like ENDLE subroutine
    end
    close(ffile)
end