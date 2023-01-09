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
        track, _ = trackfit(stubs, modules)
        # mille() per ogni hit
        for (s, m) in zip(stubs, modules)
            mille!(glder, inder, s, m, track)
        end
        # scrivi su file
        write(ffile, Int32(S+S), glder, inder) # like ENDLE subroutine
    end
    close(ffile)
end


function residuals(; ifname::String , mfname::String, title::String)
    ROOT = pyimport("ROOT")
    
    axislabels = "(LocalX_{pred} - LocalX_{hit}) [#mu m]; #Events/5 #mu m"
    
    link = ["0", "1", "2", "3", "4", "5"]
    type = ["X", "Y", "U", "V", "X", "Y"]
    res = zeros(Float32, 6)

    tree = LazyTree(ROOTFile(ifname), "Cereal", ["LocalX", "LocalY", "Link", "Bend"])
    rootfile = ROOT.TFile(replace(title, " " => "_", "-" => "_")*".root", "recreate")

    hists = []
    histchi2 = ROOT.TH1D("chi2", title*" - #chi^2; #chi^2/ndof; #Events", 200, 0, 20)
    histmedian = ROOT.TH1D("median", title * "; median" * axislabels, 400, -1000, 1000)

    for (l, t) in zip(link, type)
        histtitle = title * " - Link: " * l * " - Module type: " * t * ";"
        push!(hists, ROOT.TH1D("residuals"*l*t, histtitle*axislabels, 400, -1000, 1000))
    end
    
    modules = getmodules(mfname)
    
    stubs = StubSet{Float32}()

    for event in ProgressBar(tree)
        s = selectevent!(stubs, event)
        if !s
            continue
        end
        track, chi2 = trackfit(stubs, modules)
        histchi2.Fill(chi2/8)
        for (i, s, m, h) in zip(1:6, stubs, modules, hists)
            z = intersection(m, track)
            res[i] = -10000 * rmeas(s, z, m, track)[1] #convert in micrometers
            h.Fill(res[i])
        end
        histmedian.Fill(median(res))
    end

    histmedian.Write()
    histchi2.Write()
    for h in hists
        h.Write()
    end

    rootfile.Close()
end
