function chisquare(stubs::StubSet, modules::MUonEStation, track::Track{T}, weight::Real; cic=false) where T
    Fx = zero(T)
    Fy = zero(T)
    Fb = zero(T)
    for (s, m) in zip(stubs, modules)
        zs, zc = intersection(m, track)
        ls, lc = stub_to_local(s, m)
        rx, ry = rmeas(ls, zs, m, track)
        rb, _ = rmeas(lc, zc, m, track)
        Fx += rx^2
        Fy += ry^2
        Fb += rb^2
    end
    σx, σy, σb = sigma(weight)
    if cic
        return (Fx/σx^2 + Fy/σy^2 + Fb/σb)/2
    else
        return (Fx/σx^2 + Fb/σb)/2
    end
end

function gradient!(G, stubs::StubSet, modules::MUonEStation, track::Track{T}, weight::Real; cic=false) where T
    gx = zeros(T, 4)
    gy = zeros(T, 4)
    gb = zeros(T, 4)
    for (s, m) in zip(stubs, modules)
        zs, zc = intersection(m, track)
        ls, lc = stub_to_local(s, m)
        rx, ry = rmeas(ls, zs, m, track)
        dx, dy = derlc(zs, m)
        rb, _ = rmeas(lc, zc, m, track)
        db, _ = derlc(zc, m)

        gx -= dx * rx 
        gy -= dy * ry 
        gb -= db * rb 
    end
    σx, σy, σb = sigma(weight)
    if cic
        G .= gx/σx^2 + gy/σy^2 + gb/σb^2
    else
        G .= gx/σx^2 + gb/σb^2
    end
end

function chisquare_gradient!(F, G, stubs::StubSet, modules::MUonEStation, track::Track{T}, weight::Real; cic=false) where T
    zzs = zeros(T, 6)
    zzc = zeros(T, 6)
    rrx = zeros(T, 6)
    rry = zeros(T, 6)
    rrb = zeros(T, 6)
    vx, vy, vb = sigma(weight).^2
    for (i, s, m) in zip(1:6, stubs, modules)
        zzs[i], zzc[i] = intersection(m, track)
        ls, lc = stub_to_local(s, m)
        rrx[i], rry[i] = rmeas(ls, zzs[i], m, track)
        rrb[i], _ = rmeas(lc, zzc[i], m, track)
    end
    if G !== nothing
        gx = zeros(T, 4)
        gy = zeros(T, 4)
        gb = zeros(T, 4)
        for (zs, zc, rx, ry, rb, m) in zip(zzs, zzc, rrx, rry, rrb, modules)
            dx, dy = derlc(zs, m)
            db, _ = derlc(zc, m)

            gx -= dx * rx
            gy -= dy * ry
            gb -= db * rb
        end
        if cic
            G .= gx/vx + gy/vy + gb/vb
        else 
            G .= gx/vx + gb/vb
        end
    end
    if F !== nothing
        Fx = zero(T)
        Fy = zero(T)
        Fb = zero(T)
        for (rx, ry, rb) in zip(rrx, rry, rrb)
            Fx += rx^2
            Fy += ry^2
            Fb += rb^2
        end
        if cic
            return (Fx/vx + Fy/vy + Fb/vb)/2
        else
            return (Fx/vx + Fb/vb)/2
        end
    end        
end

function trackfit(stubs::StubSet, modules::MUonEStation, weight::Real; cic=false)
    t_0 = interpolate(stubs, modules)
    p_0 = [t_0.t0.x, t_0.t0.y, t_0.et.x, t_0.et.y]

    #f(x) = chisquare(stubs, modules, Track(x...))
    #g!(G, x) = gradient!(G, stubs, modules, Track(x...))
    fg!(F, G, x) = chisquare_gradient!(F, G, stubs, modules, Track(x...), weight, cic=cic)

    results = optimize(Optim.only_fg!(fg!), p_0, BFGS(linesearch = BackTracking()))

    popt = Optim.minimizer(results)
    #chi2 = Optim.minimum(results)
    return Track(popt...)
end

