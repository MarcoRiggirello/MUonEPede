function residuals(; nevents::Integer, mcfname::String, nmfname::String, histtitle::String)
    ROOT = pyimport("ROOT")
    histname = ("X1", "Y1", "U", "V", "X2", "Y2")

    rootfile = ROOT.TFile(histfname, "recreate")
    axislabels = ";measured hit - expected hit [cm]; #Events/0.001cm"
    hist = ROOT.TH1D("residual", histtitle*axislabels, 400, -0.2, 0.2)
    
    mcmodules = getmodules(mcfname)
    nmmodules = getmodules(nmfname)
    
    # service arrays
    stubs = StubSet{Float32}()

    for _ in ProgressBar(1:nevents)
        mcdata!(stubs, mcmodules)
        # costruisci la traccia target
        track = trackfit(stubs, nmmodules)
        # the U module
        z = intersection(nmmodules[3], track)
        res = rmeas(stubs[3], z, nmmodules[3], track)[1]
        hist.Fill(res)
    end
    hist.Write()
    rootfile.Close()
end
