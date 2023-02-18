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

function align(; ntuple::String, structure::String, iterations=4)
    tree = LazyTree(ROOTFile(ifname), "cbmsim", ["LocalX", "LocalY", "Link", "Bend"])
    
end
# align(; ntuple::String, structure::String, outcorrection?, iterations=4)
# generatebin(; ntuple::LazyTree, modules::MUonEStructure, binfile::String, weight::Union{Bool, Real})
function generatebin(; tree::LazyTree, modules::MUonEStation, binfile::String; cic=false)
    ffile = FortranFile(binfile, "w")
    # service array
    stubs = StubSet{Float32}()

    # see https://gitlab.desy.de/claus.kleinwort/millepede-ii/-/blob/main/mille.f90#L27
    # without CIC (1 rmeas + 1 sigma + 4 lder + 6 gder) * 12 measurements + 1 line of zeros
    # with CIC (1 rmeas + 1 sigma + 4 lder + 6 gder) * 18 measurements + 1 line of zeros
    S = cic ? 217 : 145
    glder = Vector{Float32}(undef, S)
    inder = Vector{Int32}(undef, S)

    glder[1] = zero(Float32)
    inder[1] = zero(Int32)
    for event in ProgressBar(tree)
        s = selectevent!(stubs, event)
        if !s
            continue
        end

        track, χ2 = trackfit(stubs, modules, cic=cic)
        ndof = cic ? 14 : 8
        weight = √(χ2 / ndof)

        for (s, m) in zip(stubs, modules)
            mille!(glder, inder, s, m, track, weight, cic=cic)
        end
        write(ffile, Int32(S + S), glder, inder) # like ENDLE subroutine
    end
    close(ffile)
end
#=
"""
    
    generatebin(; ifname::String, mfname::String, ofname::String, cfname::String, weight::Real)

Returns a Fortran binary file using as input a MUonE NTuple (vector format)
and the XML structure file.

# Arguments
- `ifname`: name of the input root file;
- `mfname`: name of the xml structure file; 
- `ofname`: name of the fortran binary output file;
- `cfname`: name of the millepede text correction file (optional);
- `weight`: weight of the u sigma;
- `cic`: if true, computes the contribution of the CIC information to the chisquare.

# Note
CAVE! The use of CIC information is, at present MUonEPede status, A WARRANTY OF BIASED ALIGNMENT.
"""
function generatebin(; ifname::String, mfname::String, ofname::String, cfnames=nothing, weight=1.0, cic=false)
    tree = LazyTree(ROOTFile(ifname), "cbmsim", ["LocalX", "LocalY", "Link", "Bend"])
    ffile = FortranFile(ofname, "w")
    modules = MUonEStation(mfname)
    if cfnames !== nothing
        for c in cfnames
            modules = applycorrections(c, modules)
        end
    end    
    # service array
    stubs = StubSet{Float32}()
    
    # see https://gitlab.desy.de/claus.kleinwort/millepede-ii/-/blob/main/mille.f90#L27
    # without CIC (1 rmeas + 1 sigma + 4 lder + 6 gder) * 12 measurements + 1 line of zeros
    # with CIC (1 rmeas + 1 sigma + 4 lder + 6 gder) * 18 measurements + 1 line of zeros
    S = cic ? 217 : 145
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
        track, chi2 = trackfit(stubs, modules, cic=cic)
        # mille() per ogni hit
        for (s, m) in zip(stubs, modules)
            mille!(glder, inder, s, m, track, weight, cic=cic)
        end
        # scrivi su file
        write(ffile, Int32(S+S), glder, inder) # like ENDLE subroutine
    end
    close(ffile)
end
=#

function residuals(; ifname::String , ofname::String, title::String, loweredges=[-200 for i in 1:6], upperedges=[200 for i in 1:6], cic=false)
    
    axislabels = "Residuals [#mu m]; Events/#mu m"
    labels2d = "Hit residuals [#mu m]; Bend residuals [#mu m]"
    
    link = ["0", "1", "2", "3", "4", "5"]
    type = ["X", "Y", "U", "V", "X", "Y"]

    s = cic ? 36 : 24
    S = cic ? 217 : 145
    uresidualindex = [i for i in 2:s:S]
    bresidualindex = [i for i in 14:s:S]
    
    binfile = FortranFile(ifname)
    rootfile = ROOT.TFile(ofname, "recreate")

    uhists = []
    bhists = []
    dhists = []

    for (l, t, le, ue) in zip(link, type, loweredges, upperedges)
        histtitle = title * "Link: " * l * " --- Type: " * t * ";"
        push!(uhists, ROOT.TH1D("residuals"*l*t, histtitle*axislabels, abs(ue - le), le, ue))
        push!(bhists, ROOT.TH1D("residuals (bend)"*l*t, histtitle*axislabels, abs(ue - le), le, ue))
        push!(dhists, ROOT.TH2D("residuals correlation"*l*t, histtitle*labels2d, abs(ue - le), le, ue, abs(ue - le), le, ue))
    end
    
    while !eof(binfile)
        _, glder, _ = read(binfile, Int32, (Float32, Int32(S)), (Int32, Int32(S)))
        for (iu, ib, uh, bh, dh) in zip(uresidualindex, bresidualindex, uhists, bhists, dhists)
		uh.Fill(10000*glder[iu])
		bh.Fill(10000*glder[ib])
        dh.Fill(10000*glder[iu], 10000*glder[ib])
        end
    end

    for h in uhists
        h.Write()
    end
    for h in bhists
        h.Write()
    end
    for h in dhists
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
