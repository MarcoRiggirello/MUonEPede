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


function applycorrections(fname::String, ms::MUonEStation)
    millepederes = readdlm(fname, skipstart=1)
    corrections = Vector{Float32}(millepederes[:,2])
    modules = Vector{MUonEModule}(undef, 6) 
    for (i, m) in enumerate(ms)
        
        x0, y0, z0 = m.r0
        θx, θy, θz = Rotations.params(m.R) 
        id = m.id
        name = m.name
        spacing = m.spacing
        j = 6 * (i - 1)
        Δx0, Δy0, Δz0 = corrections[j+1:j+3]
        Δθx, Δθy, Δθz = corrections[j+4:j+6]
        modules[i] = MUonEModule(x0 + Δx0,
                                 y0 + Δy0,
                                 z0 + Δz0,
                                 θx + Δθx,
                                 θy + Δθy,
                                 θz + Δθz,
                                 id=id,
                                 name=name,
                                 spacing=spacing)
    end
    return MUonEStation(modules...)
end
     

function selectevent!(stubs::StubSet, event)
        links = event.Link
        # filtra gli eventi con 6 hit
        if length(links) != 6
            return false
        end
        # solo 6 hit da tutti e 6 i moduli
        if sort(links) != [i for i in 0:5]
            return false
        end
        for (l, x, y, b) in zip(links, event.LocalX, event.LocalY, event.Bend)
            stubs[l+1] = Stub(x, y, b, l)
        end
        return true
end


"""
    
    generatebin(; ifname::String, mfname::String, ofname::String, cfname::String, weight::Real)

Returns a Fortran binary file using as input a MUonE NTuple (vector format)
and the XML structure file.

#Arguments
- `ifname`: name of the input root file;
- `mfname`: name of the xml structure file; 
- `ofname`: name of the fortran binary output file;
- `cfname`: name of the millepede text correction file (optional);
- `weight`: weight of the u sigma.
"""
function generatebin(; ifname::String, mfname::String, ofname::String, cfnames=nothing, weight=1.0)
    tree = LazyTree(ROOTFile(ifname), "cbmsim", ["LocalX", "LocalY", "Link", "Bend"])
    ffile = FortranFile(ofname, "w")
    modules = getmodules(mfname)
    if cfnames !== nothing
        for c in cfnames
            modules = applycorrections(c, modules)
        end
    end    
    # service array
    stubs = StubSet{Float32}()
    
    # see https://gitlab.desy.de/claus.kleinwort/millepede-ii/-/blob/main/mille.f90#L27
    S = 145 # (1 rmeas + 1 sigma + 4 lder + 6 gder) * 12 measurements + 1 line of zeros
    glder = Vector{Float32}(undef, S)
    inder = Vector{Int32}(undef, S)

    glder[1] = zero(Float32)
    inder[1] = zero(Int32)
    for event in ProgressBar(tree)
        s = selectevent!(stubs, event)
        if !s
            continue
        end
        # costruisci la traccia target
        track = trackfit(stubs, modules, weight)
        # mille() per ogni hit
        for (s, m) in zip(stubs, modules)
            mille!(glder, inder, s, m, track, weight)
        end
        # scrivi su file
        write(ffile, Int32(S+S), glder, inder) # like ENDLE subroutine
    end
    close(ffile)
end


function residuals(; ifname::String , ofname::String, title::String, loweredges=[-200 for i in 1:6], upperedges=[200 for i in 1:6])
    ROOT = pyimport("ROOT")
    
    axislabels = "Residuals [#mu m]; Events/#mu m"
    
    link = ["0", "1", "2", "3", "4", "5"]
    type = ["X", "Y", "U", "V", "X", "Y"]

    S = 145
    uresidualindex = [i for i in 2:24:S]

    binfile = FortranFile(ifname)
    rootfile = ROOT.TFile(ofname, "recreate")

    hists = []

    for (l, t, le, ue) in zip(link, type, loweredges, upperedges)
        histtitle = title * "Link: " * l * " --- Type: " * t * ";"
        push!(hists, ROOT.TH1D("residuals"*l*t, histtitle*axislabels, abs(ue - le), le, ue))
    end
    
    while !eof(binfile)
        _, glder, _ = read(binfile, Int32, (Float32, Int32(S)), (Int32, Int32(S)))
        for (i, h) in zip(uresidualindex, hists)
            h.Fill(10000*glder[i])
        end
    end

    for h in hists
        h.Write()
    end

    rootfile.Close()
end


function eigengetter(fname::String)
    f = open(fname)
    strings = readlines(f)
    eigenvalues = Vector{Float32}(undef,0)
    eigenvector = Vector{Float32}(undef,0)
    eigenvectors = []
    for s in strings
        ss = split(s, r" +", keepempty=false)
        ll = length(ss)       
        if ll == 5
            append!(eigenvalues, parse(Float32, ss[end]))
        elseif ll == 0
            append!(eigenvectors, eigenvector)
            eigenvector = []
        elseif  ll == 6 || ll == 4 || ll == 2
            for i in eachindex(ss)
                if iseven(i)
                    append!(eigenvector, parse(Float32, ss[i]))
                end
            end
        else
            throw("not expected.")
        end
    end
    return Diagonal(eigenvalues), hcat(eigenvectors...)
end
