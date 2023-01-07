function residualsmc(; nevents::Integer, mcfname::String, nmfname::String, title::String)
    ROOT = pyimport("ROOT")
    
    axislabels = ";(measured LocalX - fitted LocalX) [cm]; #Events/0.001cm"
    
    link = ["0", "1", "2", "3", "4", "5"]
    type = ["X", "Y", "U", "V", "X", "Y"]

    hists = []
    rootfile = ROOT.TFile(replace(title, " " => "_", "-" => "_")*".root", "recreate")

    for (l, t) in zip(link, type)
        histtitle = title * " - Link: " * l * " - Module type: " * t * " (MC)"
        push!(hists, ROOT.TH1D("residuals"*l*t, histtitle*axislabels, 400, -0.2, 0.2))
    end
    
    mcmodules = getmodules(mcfname)
    nmmodules = getmodules(nmfname)
    
    stubs = StubSet{Float32}()

    for _ in ProgressBar(1:nevents)
        mcdata!(stubs, mcmodules)
        track = trackfit(stubs, nmmodules)
        for (s, m, h) in zip(stubs, nmmodules, hists)
            z = intersection(m, track)
            res = rmeas(s, z, m, track)[1]
            h.Fill(res)
        end
    end

    for h in hists
        h.Write()
    end

    rootfile.Close()
end
