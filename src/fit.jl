function chisquare_gradient!(F, G, stubs::StubSet, modules::MUonEStation, track::Track{T}; cic::Bool, skip::Union{Bool, Integer}) where T
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
        if s.link == skip
            continue
        end
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

# approximation!
function covariance_matrix!(H, stubs::StubSet, modules::MUonEStation, track::Track{T}; cic::Bool, skip::Union{Bool, Integer}) where T
    H .= 0
    for (s, m) in zip(stubs, modules)
        if s.link == skip
            continue
        end
        zs, zc = intersection(m, track)
        σ2x, σ2y, σ2b = sigma(s).^2
        dx, dy = derlc(zs, m)
        db, _  = derlc(zc, m)
        for i in 1:4
            for j in 1:4
                H[i,j] += dx[i] * dx[j] / σ2x
                H[i,j] += db[i] * db[j] / σ2b
                if cic
                    H[i,j] += dy[i] * dy[j] / σ2y
                end
            end
        end
    end
    H .= inv(H)
end


function trackfit(stubs::StubSet, modules::MUonEStation; cic=false, skip=false)
    t_0 = interpolate(stubs, modules)
    p_0 = [t_0.t0.x, t_0.t0.y, t_0.et.x, t_0.et.y]

    fg!(F, G, x) = chisquare_gradient!(F, G, stubs, modules, Track(x...), cic=cic, skip=skip)
    
    local results

    n_fit_max = 50

    for i in 0:n_fit_max
        results = optimize(Optim.only_fg!(fg!),
                           p_0,
                           LBFGS(P = diagm([1f-2, 1f-2, 1f2, 1f2]),
                                 linesearch = LineSearches.MoreThuente(gtol = 0.1)),
                           Optim.Options(g_tol = (i == 0) ? 1e4 : 1e-7))
        if Optim.g_converged(results) && i != 0
            break
        end
        p_0 = Optim.minimizer(results) + randn(4) .* [2e-2, 2e-2, 2e-4, 2e-4] 
        if i == n_fit_max
            @warn "gradient norm is greater than 1e-7"
        end
    end

    if Optim.converged(results)
        return results
    else
        throw("local fit not converged!")
    end
end


function trackfit!(popt, stubs::StubSet, modules::MUonEStation; cic=false, skip=false)
    results = trackfit(stubs, modules; cic=cic, skip=skip)
    popt .= Optim.minimizer(results)
    chi2 = Optim.minimum(results)
    return chi2
end


function trackfit!(popt, covm, stubs::StubSet, modules::MUonEStation; cic=false, skip=false)
    results = trackfit(stubs, modules; cic=cic, skip=skip)
    popt .= Optim.minimizer(results)
    covariance_matrix!(covm, stubs, modules, Track(popt...), cic=cic, skip=skip)
    chi2 = Optim.minimum(results)
    return chi2
end

