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
    return MUonEStation(modules...)
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
    tree = LazyTree(ROOTFile(ifname), "Cereal", ["LocalX", "LocalY", "Link", "Bend"])
    ffile = FortranFile(ofname, "w")
    modules = getmodules(mfname)
    
    # service array
    stubs = Vector{Stub{Float32}}(undef, 6)

    # see https://gitlab.desy.de/claus.kleinwort/millepede-ii/-/blob/main/mille.f90#L27
    S = 145 # (1 rmeas + 1 sigma + 4 lder + 6 gder) * 12 measurements + 1 line of zeros
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
        # solo 6 hit da tutti e 6 i moduli
        if sort(links) != [i for i in 0:5]
            continue
        end
        for (l, x, y, b) in zip(links, event.LocalX, event.LocalY, event.Bend)
            stubs[l+1] = Stub(x, y, b, l)
        end
        # costruisci la traccia target
        track = trackfit(stubs, modules)
        # mille() per ogni hit
        for (s, m) in zip(stubs, modules)
            mille!(glder, inder, s, m, track)
        end
        # scrivi su file
        write(ffile, Int32(S+S), glder, inder) # like ENDLE subroutine
    end
    close(ffile)
end