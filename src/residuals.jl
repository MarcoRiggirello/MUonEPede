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
