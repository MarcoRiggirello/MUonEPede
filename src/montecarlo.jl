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
        track, _ = trackfit(stubs, nmmodules)
        # mille() per ogni hit
        for (s, m) in zip(stubs, nmmodules)
            mille!(glder, inder, s, m, track)
        end
        # scrivi su file
        write(ffile, Int32(S+S), glder, inder) # like ENDLE subroutine
    end
    close(ffile)
end


function residualsmc(; nevents::Integer, mcfname::String, nmfname::String, title::String)
    ROOT = pyimport("ROOT")
    
    axislabels = "(LocalX#_{pred} - LocalX#_{hit}) [#mu m]; #Events/5 #mu m"
    
    link = ["0", "1", "2", "3", "4", "5"]
    type = ["X", "Y", "U", "V", "X", "Y"]
    res = zeros(Float32, 6)

    rootfile = ROOT.TFile(replace(title, " " => "_", "-" => "_")*".root", "recreate")

    hists = []
    histchi2 = ROOT.TH1D("chi2", title*" - #chi#^2; #chi#^2; #Events", 200, 0, 100)
    histmedian = ROOT.TH1D("median", title * "; median" * axislabels, 400, -1000, 1000)

    for (l, t) in zip(link, type)
        histtitle = title * " - Link: " * l * " - Module type: " * t * ";"
        push!(hists, ROOT.TH1D("residuals"*l*t, histtitle*axislabels, 400, -1000, 1000))
    end
    
    mcmodules = getmodules(mcfname)
    nmmodules = getmodules(nmfname)
    
    stubs = StubSet{Float32}()

    for _ in ProgressBar(1:nevents)
        mcdata!(stubs, mcmodules)
        track, chi2 = trackfit(stubs, nmmodules)
        histchi2.Fill(chi2)
        for (i, s, m, h) in zip(1:6, stubs, nmmodules, hists)
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
