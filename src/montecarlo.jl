# a dummy Monte Carlo for pede alignment

"""
   Do the inverse of strip_to_local 
"""
function local_to_stub(q::StaticVector{3,T}, m::MUonEModule) where {T<:Real}
    nstrips = 1016
    strip_pitch = 0.009

    strip_X = 2 * q.x ÷ strip_pitch + nstrips + 3
    strip_Y = q.y > 0 ? 0.75 : 0.25
    # we don't use the bend at the moment
    return Stub{T}(strip_X, strip_Y, 0.0, m.id)
end

function mcdata!(stubs::StubSet, modules::MUonEStation)
    x0 = randn(Float32)
    y0 = randn(Float32)
    mx = 1f-3 * randn(Float32)
    my = 1f-3 * randn(Float32)
    
    t = Track{Float32}(x0, y0, mx, my)
    z = intersection.(modules, [t])
    
    stubs[1:6] = local_to_stub.(global_to_local.(t.(z), modules), modules)
end

"""
    generatebinmc(; nevents::Integer, mcfname::String, nmfname::String, ofname::String)

Generate MonteCarlo data in a Fortran binary file using 
as input the XML structure file.

#Arguments
- `nevents`: number of tracks;
- `mcfname`: name of the xml structure file with the montecarlo truth;
- `nmfnmae`: name of the xml structure file with nominal positions; 
- `ofname`: name of the fortran binary output file.
"""
function generatebinmc(; nevents::Integer, mcfname::String, nmfname::String, ofname::String)
    ffile = FortranFile(ofname, "w")
    mcmodules = getmodules(mcfname)
    nmmodules = getmodules(nmfname)
    
    # service array
    stubs = StubSet{Float32}()

    # see https://gitlab.desy.de/claus.kleinwort/millepede-ii/-/blob/main/mille.f90#L27
    S = 145 # (1 rmeas + 1 sigma + 4 lder + 6 gder) * 12 measurements + 1 line of zeros
    glder = Vector{Float32}(undef, S)
    inder = Vector{Int32}(undef, S)

    glder[1] = zero(Float32)
    inder[1] = zero(Int32)
    for _ in ProgressBar(1:nevents)
        mcdata!(stubs, mcmodules)
        # costruisci la traccia target
        track = trackfit(stubs, nmmodules)
        # mille() per ogni hit
        for (s, m) in zip(stubs, nmmodules)
            mille!(glder, inder, s, m, track)
        end
        # scrivi su file
        write(ffile, Int32(S+S), glder, inder) # like ENDLE subroutine
    end
    close(ffile)
end