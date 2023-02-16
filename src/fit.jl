function chisquare(stubs::StubSet, modules::MUonEStation, track::Track{T}; cic=false) where T
    Fx = zero(T)
    Fy = zero(T)
    Fb = zero(T)
    for (s, m) in zip(stubs, modules)
        zs, zc = intersection(m, track)
        ls, lc = stub_to_local(s, m)
        rx, ry = rmeas(ls, zs, m, track)
        rb, _  = rmeas(lc, zc, m, track)
        σx, σy, σb = sigma(s)
        Fx += (rx/σx)^2
        Fy += (ry/σy)^2
        Fb += (rb/σb)^2
    end
    if cic
        return (Fx + Fy + Fb)/2
    else
        return (Fx + Fb)/2
    end
end


function chisquare_gradient!(F, G, stubs::StubSet, modules::MUonEStation, track::Track{T}; cic=false) where T
    f = F !== nothing
    g = G !== nothing

    if g
        Gx = zeros(T, 4)
        Gy = zeros(T, 4)
        Gb = zeros(T, 4)
    end

    if f
        Fx = zero(T)
        Fy = zero(T)
        Fb = zero(T)
    end

    for (s, m) in zip(stubs, modules)
        zs, zc = intersection(m, track)
        ls, lc = stub_to_local(s, m)

        rx, ry = rmeas(ls, zs, m, track)
        rb, _  = rmeas(lc, zc, m, track)

        σ2x, σ2y, σ2b = sigma(s).^2

        if g
            dx, dy = derlc(zs, m)
            db, _  = derlc(zc, m)

            Gx -= dx * rx / σ2x
            Gy -= dy * ry / σ2y
            Gb -= db * rb / σ2b
        end

        if f
            Fx += rx^2 / σ2x
            Fy += ry^2 / σ2y
            Fb += rb^2 / σ2b
        end
    end
 
    if g
        G .= cic ? Gx + Gy + Gb : Gx + Gb 
    end
   
    if f
        return cic ? (Fx + Fy + Fb)/2 : (Fx + Fb)/2
    end        
end

function trackfit(stubs::StubSet, modules::MUonEStation; cic=false)
    #t_0 = interpolate(stubs, modules)
    #p_0 = [t_0.t0.x, t_0.t0.y, t_0.et.x, t_0.et.y]
    p_0 = [0.0 for i in 1:4]

    #f = TwiceDifferentiable(x->chisquare(stubs, modules, Track(x...), weight), p_0, autodiff=:forward)
    fg!(F, G, x) = chisquare_gradient!(F, G, stubs, modules, Track(x...), cic=cic)

    results = optimize(Optim.only_fg!(fg!),
                       p_0,
                       BFGS(alphaguess = InitialHagerZhang()))
    popt = Optim.minimizer(results)
    chi2 = Optim.minimum(results)
    if Optim.converged(results)
        return Track(popt...), chi2
    else
        throw("local fit not converged!")
    end
end

