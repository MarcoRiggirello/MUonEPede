function chisquare(stubs::StubSet, modules::MUonEStation, track::Track{T}) where T
    Fx = zero(T)
    Fy = zero(T)
    for (s, m) in zip(stubs, modules)
        z = intersection(m, track)
        rx, ry = rmeas(s, z, m, track)
        Fx += rx^2
        Fy += ry^2
    end
    sx, sy = sigma()
    return (Fx/sx^2 + Fy/sy^2)/2
end

function gradient!(G, stubs::StubSet, modules::MUonEStation, track::Track{T}) where T
    gx = zeros(T, 4)
    gy = zeros(T, 4)
    for (s, m) in zip(stubs, modules)
        z = intersection(m, track)
        rx, ry = rmeas(s, z, m, track)
        dx, dy = derlc(z, m)

        gx -= dx * rx 
        gy -= dy * ry 
    end
    sx, sy = sigma()
    G .= gx/sx^2 + gy/sy^2
end


function trackfit(stubs::StubSet, modules::MUonEStation)
    t_0 = interpolate(stubs, modules)
    p_0 = [t_0.t0.x, t_0.t0.y, t_0.et.x, t_0.et.y]

    f(x) = chisquare(stubs, modules, Track(x...))
    g!(G, x) = gradient!(G, stubs, modules, Track(x...))

    popt = Optim.minimizer(optimize(f,
                                    g!,
                                    p_0,
                                    BFGS(linesearch = LineSearches.BackTracking())))

    return Track(popt...)
end

