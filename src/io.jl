"""
    getmodules(fname::String)
    
Returns a `SVector` of `MUonEModules` constructed from an xml file.
See https://gitlab.cern.ch/muesli/daq-sw/daq-decode/-/blob/stdVec_bxAssembled/Structure/MUonEStructure_TB2022.xml
"""
function getmodules(fname::String)
    doc = readxml(fname)
    station = doc.root.firstelement
    modulenodes = findall("//Module", station)
    modules = Vector{MUonEModule}(undef, 6) 
    for (i, m) in enumerate(modulenodes)
        pos, rot = elements(m)
        
        x0 = parse(Float32, pos["positionX"])
        y0 = parse(Float32, pos["positionY"])
        z0 = parse(Float32, pos["positionZ"])
        θx = parse(Float32, rot["rotationX"])
        θy = parse(Float32, rot["rotationY"])
        θz = parse(Float32, rot["rotationZ"])
        id = parse(Int32, m["linkId"])
        name = m["name"]
        spacing = parse(Float32, m["sensorSpacing"])

        modules[i] = MUonEModule(x0, y0, z0, θx, θy, θz, id=id, name=name, spacing=spacing)
    end
    return SVector{6}(modules)
end


"""
    generatebin(; ifname::String, mfname::String, ofname::String)

Returns a Fortran binary file using as input a MUonE NTuple (vector format)
and the XML structure file.

#Arguments
- `ifname::String`: name of the input root file;
- `mfname::String`: name of the xml structure file; 
- `ofname::String`: name of the fortran binary output file.
"""
function generatebin(; ifname::String, mfname::String, ofname::String)
    tree = LazyTree(ROOTFile(ifname), "Cereal", ["Link", "LocalX", "LocalY"])
    ffile = FortranFile(ofname, "w")
    modules = getmodules(mfname)
    
    # service arrays
    ix = MVector{6, Int32}(undef)
    strips_X = MVector{6, Float32}(undef)
    strips_Y = MVector{6, Float32}(undef)

    # see https://gitlab.desy.de/claus.kleinwort/millepede-ii/-/blob/main/mille.f90#L27
    S = 73 # (1 rmeas + 1 sigma + 4 lder + 6 gder) * 6 measurements + 1 line of zeros
    glder = Vector{Float32}(undef, S)
    inder = Vector{Int32}(undef, S)

    glder[1] = zero(Float32)
    inder[1] = zero(Int32)
    for event in ProgressBar(tree)
        links = event.Link
        # filtra gli eventi con 6 hit
        if length(links) != 6
            continue
        end
        sortperm!(ix, links)
        # solo 6 hit da tutti e 6 i moduli
        if links[ix] != [i for i in 0:5]
            continue
        end
        strips_X[1:6] = event.LocalX
        strips_Y[1:6] = event.LocalY
        permute!(strips_X, ix)
        permute!(strips_Y, ix)
        # costruisci la traccia target
        hit0 = local_to_global(strip_to_local(strips_X[1], strips_Y[1], modules[1]), modules[1])
        hit1 = local_to_global(strip_to_local(strips_X[2], strips_Y[2], modules[2]), modules[2])
        hit4 = local_to_global(strip_to_local(strips_X[5], strips_Y[5], modules[5]), modules[5])
        hit5 = local_to_global(strip_to_local(strips_X[6], strips_Y[6], modules[6]), modules[6])
        track = interpolate(hit0, hit1, hit4, hit5)
        # mille() per ogni hit
        for (sx, sy, m) in zip(strips_X, strips_Y, modules)
            mille!(glder, inder, sx, sy, m, track)
        end
        # scrivi su file
        write(ffile, Int32(S+S), glder, inder) # like ENDLE subroutine
    end
    close(ffile)
end