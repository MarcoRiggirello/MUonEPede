function residuals(; nevents::Integer, mcfname::String, nmfname::String, histfname::String, histtitle="residuals")
    ROOT = pyimport("ROOT")

    rootfile = ROOT.TFile(histfname, "recreate")
    axislabels = ";measured hit - expected hit [cm]; #Events/0.001cm"
    hist = ROOT.TH1D("residual", histtitle*axislabels, 400, -0.2, 0.2)
    
    mcmodules = getmodules(mcfname)
    nmmodules = getmodules(nmfname)
    
    # service arrays
    stubs = Vector{Stub{Float32}}(undef, 6)

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
