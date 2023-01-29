# a dummy Monte Carlo for pede alignment

"""
   Do the inverse of strip_to_local 
"""
function local_to_stub(q::StaticVector{3,T}, m::MUonEModule) where {T<:Real}
    nstrips = 1016
    strip_pitch = 0.009

    strip_X = q.x / strip_pitch + nstrips/2 - 1/2
    strip_Y = q.y > 0 ? 0.75 : 0.25
    # we don't use the bend at the moment
    return Stub{T}(strip_X, strip_Y, 0.0, m.id)
end

function mcdata!(stubs::StubSet, modules::MUonEStation; beamsize_x=1.0, beamsize_y=1.0)
    while true
        x0 = Float32(beamsize_x) * randn(Float32)
        y0 = Float32(beamsize_y) * randn(Float32)
        mx = 1f-3 * randn(Float32)
        my = 1f-3 * randn(Float32)
    
        t = Track{Float32}(x0, y0, mx, my)
        z = intersection.(modules, [t])

        qq = global_to_local.(t.(z), modules)
        if all(q -> abs(q[1] < 5), qq) && all(q -> abs(q[2] < 5), qq)
            stubs[1:6] = local_to_stub.(qq, modules)
            break
        end
    end
end

"""
    generatebinmc(; nevents::Integer, mcfname::String, nmfname::String, ofname::String, cfname=nothing, weight=1.0, beamsize_x=1.0, beamsize_y=1.0)

Generate MonteCarlo data in a Fortran binary file using 
as input the XML structure file.

#Arguments
- `nevents`: number of tracks;
- `mcfname`: name of the xml structure file with the montecarlo truth;
- `nmfnmae`: name of the xml structure file with nominal positions; 
- `ofname`: name of the fortran binary output file;
- `cfnames`: list of parameters correction (the expected syntax is the one of the typical millepede.res file);
- `weight`: weight of sigma;
- `cic`: use CIC information in chisquare;
- `beamsize_x`: the x sigma of the beam (in cm);
- `beamsize_y`: the y sigma of the beam (in cm).
"""
function generatebinmc(; nevents::Integer, mcfname::String, nmfname::String, ofname::String, cfnames=nothing, weight=1.0, cic=false, beamsize_x=1.0, beamsize_y=1.0)
    ffile = FortranFile(ofname, "w")
    mcmodules = getmodules(mcfname)
    nmmodules = getmodules(nmfname)
    if cfnames !== nothing
        for c in cfnames
            nmmodules = applycorrections(c, nmmodules)
        end
    end    
    
    # service array
    stubs = StubSet{Float32}()

    # see https://gitlab.desy.de/claus.kleinwort/millepede-ii/-/blob/main/mille.f90#L27
    # with CIC: (1 rmeas + 1 sigma + 4 lder + 6 gder) * 12 measurements + 1 line of zeros
    # without CIC: (1 rmeas + 1 sigma + 4 lder + 6 gder) * 6 measurements + 1 line of zeros
    S = cic ? 145 : 73
    glder = Vector{Float32}(undef, S)
    inder = Vector{Int32}(undef, S)

    glder[1] = zero(Float32)
    inder[1] = zero(Int32)
    for _ in ProgressBar(1:nevents)
        mcdata!(stubs, mcmodules, beamsize_x=beamsize_x, beamsize_y=beamsize_y)
        # costruisci la traccia target
        track = trackfit(stubs, nmmodules, weight, cic=cic)
        # mille() per ogni hit
        for (s, m) in zip(stubs, nmmodules)
            mille!(glder, inder, s, m, track, weight, cic=cic)
        end
        # scrivi su file
        write(ffile, Int32(S+S), glder, inder) # like ENDLE subroutine
    end
    close(ffile)
end
