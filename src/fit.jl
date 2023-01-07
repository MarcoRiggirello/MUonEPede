function chisquare(stubs::StubSet, modules::MUonEStation, track::Track{T}) where T
    F = zero(T)
    for (s, m) in zip(stubs, modules)
        z = intersection(m, track)
        rx, ry = rmeas(s, z, m, track)
        sx, sy = sigma()
        F += (rx/sx)^2 + (ry/sy)^2
    end
    return F/2
end

function gradient!(G, modules::MUonEStation, track::Track)
    G = zero(G)
    for m in modules
        z = intersection(m, track)
        dx, dy = derlc(z, m)
        G -= dx + dy
    end
end

function trackfit(stubs::StubSet, modules::MUonEStation)
    t_0 = interpolate(stubs, modules)
    p_0 = [t_0.t0.x, t_0.t0.y, t_0.et.x, t_0.et.y]

    f(x) = chisquare(stubs, modules, Track(x...))
    g!(G, x) = gradient!(G, modules, Track(x...))

    popt = Optim.minimizer(optimize(f, g!, p_0, BFGS()))

    return Track(popt...)
end

