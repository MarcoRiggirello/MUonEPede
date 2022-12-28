function residuals(; nevents::Integer, mcfname::String, nmfname::String, histfname::String, histtitle="residuals")
    
    rootfile = ROOT.TFile(histfname, "recreate")
    axislabels = ";measured hit - expected hit [cm]; #Events/0.01cm"
    hist = ROOT.TH1D("residual", histtitle*axislabels, 400, -2., 2.)
    
    mcmodules = getmodules(mcfname)
    nmmodules = getmodules(nmfname)
    
    # service arrays
    strips_X = MVector{6, Float32}(undef)
    strips_Y = MVector{6, Float32}(undef)

    for _ in ProgressBar(1:nevents)
        mcdata!(strips_X, strips_Y, mcmodules)
        # costruisci la traccia target
        hit0 = local_to_global(strip_to_local(strips_X[1], strips_Y[1], nmmodules[1]), nmmodules[1])
        hit1 = local_to_global(strip_to_local(strips_X[2], strips_Y[2], nmmodules[2]), nmmodules[2])
        hit4 = local_to_global(strip_to_local(strips_X[5], strips_Y[5], nmmodules[5]), nmmodules[5])
        hit5 = local_to_global(strip_to_local(strips_X[6], strips_Y[6], nmmodules[6]), nmmodules[6])
        track = interpolate(hit0, hit1, hit4, hit5)

        # the U module
        z = intersection(nmmodules[3], track)
        res = rmeas(strips_X[3], strips_Y[3], z, nmmodules[3], track)
        hist.Fill(res)
    end
    hist.Write()
    rootfile.Close()
end
