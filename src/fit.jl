function chisquare(stubs::StubSet, modules::MUonEStation, track::Track{T}, weight::Real) where T
    Fx = zero(T)
    Fy = zero(T)
    for (s, m) in zip(stubs, modules)
        z = intersection(m, track)
        rx, ry = rmeas(s, z, m, track)
        Fx += rx^2
        Fy += ry^2
    end
    sx, sy = sigma(weight)
    return (Fx/sx^2 + Fy/sy^2)/2
end

function gradient!(G, stubs::StubSet, modules::MUonEStation, track::Track{T}, weight::Real) where T
    gx = zeros(T, 4)
    gy = zeros(T, 4)
    for (s, m) in zip(stubs, modules)
        z = intersection(m, track)
        rx, ry = rmeas(s, z, m, track)
        dx, dy = derlc(z, m)

        gx -= dx * rx 
        gy -= dy * ry 
    end
    sx, sy = sigma(weight)
    G .= gx/sx^2 + gy/sy^2
end

function chisquare_gradient!(F, G, stubs::StubSet, modules::MUonEStation, track::Track{T}, weight::Real) where T
    zz = zeros(T, 6)
    rrx = zeros(T, 6)
    rry = zeros(T, 6)
    vx, vy = sigma(weight).^2
    for (i, s, m) in zip(1:6, stubs, modules)
        zz[i] = intersection(m, track)
        rrx[i], rry[i] = rmeas(s, zz[i], m, track)
    end
    if G !== nothing
        gx = zeros(T, 4)
        gy = zeros(T, 4)
        for (z, rx, ry, m) in zip(zz, rrx, rry, modules)
            dx, dy = derlc(z, m)

            gx -= dx * rx
            gy -= dy * ry
        end
        G .= gx/vx + gy/vy
    end
    if F !== nothing
        Fx = zero(T)
        Fy = zero(T)
        for (rx, ry) in zip(rrx, rry)
            Fx += rx^2
            Fy += ry^2
        end
        return (Fx/vx + Fy/vy)/2
    end        
end

function trackfit(stubs::StubSet, modules::MUonEStation, weight::Real)
    t_0 = interpolate(stubs, modules)
    p_0 = [t_0.t0.x, t_0.t0.y, t_0.et.x, t_0.et.y]

    #f(x) = chisquare(stubs, modules, Track(x...))
    #g!(G, x) = gradient!(G, stubs, modules, Track(x...))
    fg!(F, G, x) = chisquare_gradient!(F, G, stubs, modules, Track(x...), weight)

    results = optimize(Optim.only_fg!(fg!), p_0, BFGS(linesearch = BackTracking()))

    popt = Optim.minimizer(results)
    #chi2 = Optim.minimum(results)
    return Track(popt...)
end

