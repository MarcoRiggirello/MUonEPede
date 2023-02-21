# a dummy Monte Carlo for pede alignment

"""

Do the inverse of strip_to_local 

"""
function local_to_stub(q_s::StaticVector{3,T}, q_c::StaticVector{3,T}, m::MUonEModule, offset::Real) where {T<:Real}
    nstrips = 1016
    strip_pitch = 0.009

    qsx = q_s.x + T(3.f-3) * randn(T)
    qcx = q_c.x + T(3.f-3) * randn(T)

    strip_X = round(qsx / strip_pitch, digits = 1, base = 2) + nstrips/2 - 1/2
    strip_Y = q_s.y > 0 ? 0.75 : 0.25
    strip_c = round((qcx + T(offset)) / strip_pitch, digits = 1, base = 2) + nstrips/2 - 1/2
    bend = (strip_c - strip_X)
    return Stub{T}(strip_X, strip_Y, bend, m.id)
end

function mcdata!(stubs::StubSet, modules::MUonEStation{Float32}; beamsigma_x=1.0, beamsigma_y=1.0, offsets=[0.0f0 for i in 1:6])
    while true
        x0 = Float32(beamsigma_x) * randn(Float32)
        y0 = Float32(beamsigma_y) * randn(Float32)
        mx = 1f-3 * randn(Float32)
        my = 1f-3 * randn(Float32)
    
        t = Track{Float32}(x0, y0, mx, my)
        z = intersection.(modules, [t])
        z_s = getindex.(z, 1)
        z_c = getindex.(z, 2)

        qq_s = global_to_local.(t.(z_s), modules)
        qq_c = global_to_local.(t.(z_c), modules)
        if all(q -> abs(q[1]) < 5, qq_s) && all(q -> abs(q[2]) < 5, qq_s)
            stubs[1:6] = local_to_stub.(qq_s, qq_c, modules, offsets)
            return x0, y0, mx, my
        end
    end
end


"""
    
    toymontecarlo(N::Integer; station::String, out::String, beamsigma_x=1.0, beamsigma_y=1.0, offsets=[0.0f0 for i in 1:6])

Generate MonteCarlo data using the XML structure file as input.

#Arguments
- `N`: number of tracks;
- `station`: name of the xml structure file with the montecarlo true postitions;
- `nmfnmae`: name of the xml structure file with nominal positions; 
- `out`: name of the root output file;
- `beamsize_x`: the x sigma of the beam (in cm);
- `beamsize_y`: the y sigma of the beam (in cm);
- `offsets`: misalignment for seed and correlation layer.

"""
function toymontecarlo(N::Integer; station::String, out::String, beamsigma_x=1.0, beamsigma_y=1.0, offsets=[0.0f0 for i in 1:6])
    f = ROOT.TFile(out, "recreate")
    m = MUonEStation{Float32}(station)

    t = ROOT.TTree("cbmsim", "cbmsim")

    nstubs = [0x6]
    link = [UInt16(i) for i in 0:5] 

    z = [Float32(mm.r0.z) for mm in m]
    localx = Vector{Float32}(undef,6)
    localy = Vector{Float32}(undef,6)
    bend = Vector{Float32}(undef,6)

    t.Branch("nStubs", nstubs, "nStubs/b")
    t.Branch("Link", link, "Link[nStubs]/s")
    t.Branch("Z", z, "Z[nStubs]/F")
    t.Branch("LocalX", localx, "LocalX[nStubs]/F")
    t.Branch("LocalY", localy, "LocalY[nStubs]/F")
    t.Branch("Bend", bend, "Bend[nStubs]/F")

    bx = [0x0000]
    superid = [0x00000000]

    t.Branch("Bx", bx, "Bx/s")
    t.Branch("SuperID", superid, "SuperID/i")

    mx = [0.f0]
    my = [0.f0]
    x0 = [0.f0]
    y0 = [0.f0]

    t.Branch("Track_mx", mx, "mx/F")
    t.Branch("Track_my", my, "my/F")
    t.Branch("Track_x0", x0, "x0/F")
    t.Branch("Track_y0", y0, "y0/F")

    s = StubSet{Float32}()
    
    for i in ProgressBar(1:N)
        superid[1], bx[1] = divrem(i, 3186)

        x0[1], y0[1] ,mx[1], my[1] = mcdata!(s, m, beamsigma_x=beamsigma_x, beamsigma_y=beamsigma_y, offsets=offsets)
        for (j, stub) in enumerate(s)
            localx[j] = stub.localX                     
            localy[j] = stub.localY
            bend[j] = stub.bend
        end

        t.Fill()
    end
    t.Write()
    f.Close()
end
