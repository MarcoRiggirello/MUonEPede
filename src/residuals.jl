function chisquare(stubs, modules, track::Track{T}) where T
    F = zero(T)
    for (s, m) in zip(stubs, modules)
        z = intersection(m, track)
        rx, ry = rmeas(s, z, m, track)
        sx, sy = sigma()
        F += (rx/sx)^2 + (ry/sy)^2
    end
    return F/2
end

function gradient(modules, track::Track{T}) where T
    G = zeros(T, 4)
    for m in modules
        z = intersection(m, track)
        dx, dy = derlc(z, m)
        G -= dx + dy
    end
    return G
end

function fit(stubs::AbstractVector{Stub}, modules::AbstractVector{MUonEModule})
    t_0 = interpolate(stubs, modules)
    p_0 = [t_0.t0.x, t_0.t0.y, t_0.et.x, t_0.et.y]

    f(x) = chisquare(stubs, modules, Track(x...))
    g(x) = gradient(modules, Track(x...))

    res = optimize(f, g, p_0, inplace = false)
    popt = Optim.minimizer(res)

    return Track(popt...)
end

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
        track = fit(stubs, nmmodules)
        # the U module
        z = intersection(nmmodules[3], track)
        res = rmeas(stubs[3], z, nmmodules[3], track)[1]
        hist.Fill(res)
    end
    hist.Write()
    rootfile.Close()
end
