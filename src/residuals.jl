# residui bias e non bias
# i bias: resX vs bend
# MC track spot and redisuals
# de poi atutot
function fastresiduals(; ifname::String , ofname::String, title::String, loweredges=[-200 for i in 1:6], upperedges=[200 for i in 1:6], cic=false)
    
    axislabels = "hit_meas - hit_fit [#mu m]; Events/#mu m"
    labels2d = "seed layer residuals [#mu m]; correlation layer residuals [#mu m]"
    
    link = ["0", "1", "2", "3", "4", "5"]
    type = ["X", "Y", "U", "V", "X", "Y"]

    s = cic ? 36 : 24
    S = cic ? 217 : 145
    uresidualindex = [i for i in 2:s:S]
    bresidualindex = [i for i in 14:s:S]
    
    binfile = FortranFile(ifname)
    rootfile = ROOT.TFile(ofname, "recreate")

    uhists = []
    bhists = []
    dhists = []

    for (l, t, le, ue) in zip(link, type, loweredges, upperedges)
        histtitle = title * "Link: " * l * " --- Type: " * t * ";"
        push!(uhists, ROOT.TH1D("residuals"*l*t, histtitle*axislabels, abs(ue - le), le, ue))
        push!(bhists, ROOT.TH1D("residuals (bend)"*l*t, histtitle*axislabels, abs(ue - le), le, ue))
        push!(dhists, ROOT.TH2D("residuals correlation"*l*t, histtitle*labels2d, abs(ue - le), le, ue, abs(ue - le), le, ue))
    end
    while !eof(binfile)
        _, glder, _ = read(binfile, Int32, (Float32, Int32(S)), (Int32, Int32(S)))
        for (iu, ib, uh, bh, dh) in zip(uresidualindex, bresidualindex, uhists, bhists, dhists)
		uh.Fill(10000*glder[iu])
        bh.Fill(10000*glder[ib])
        dh.Fill(10000*glder[iu], 10000*glder[ib])
        end
    end

    for h in uhists
        h.Write()
    end
    for h in bhists
        h.Write()
    end
    for h in dhists
        h.Write()
    end

    rootfile.Close()
end

function residual(s::Stub, m::MUonEModule, t::Track)
    z_s, z_c = intersection(m, t)
    l_s, l_c = stub_to_local(s, m)

    r_s, _ = rmeas(l_s, z_s, m, t)
    r_c, _ = rmeas(l_c, z_c, m, t)

    return 10_000 * r_s, 10_000 * r_c
end


function pull(s::Stub, m::MUonEModule, t::Track, w::Real, pcov::AbstractMatrix)
    zs, zc = intersection(m, t)
    ls, lc = stub_to_local(s, m)

    rs, _ = rmeas(ls, zs, m, t)
    rc, _ = rmeas(lc, zc, m, t)
    
    vs, _, vc = sigma(s, w).^2
    
    js = derlc(zs, m)
    jc = derlc(zc, m)

    σs = √(vs + dot(js, pcov, js)) # o ci va il meno??
    σc = √(vc + dot(jc, pcov, jc))
    return rs/σs, rc/σc 
end 


function residuals(; tree::LazyTree, modules::MUonEStation, rootfile::String, title::String, cic=false)
    ismontecarlo = any(i -> i == :Track_mx, propertynames(tree))

    ll_residuals = "; hit_meas - hit_fit [#mu m]; Events/#mu m"
    ll_pulls = "; (hit_meas - hit_fit)/sigma; Events"
    ll_seed_vs_corr = "; seed layer residuals [#mu m]; correlation layer residuals [#mu m]"
    ll_res_vs_bend = "; hit_meas - hit_fit [#mu m]; Bend [strips]"
    
    link = ["0", "1", "2", "3", "4", "5"]
    #type = ["X", "Y", "U", "V", "X", "Y"]

    hh_seed = [] 
    hh_corr = [] 
    hh_seed_vs_corr = [] 
    hh_res_vs_bend = []

    hh_seed_unbiased = []
    hh_corr_unbiased = []
    
    f = ROOT.TFile(rootfile, "recreate")

    e = 200
    e

    for l in link
        histtitle = title * "Link: " * l
        push!(hh_seed, ROOT.TH1D("seed" * l , histtitle * " -- Seed Layer "  * ll_residuals, 2 * e, -e, e))
        push!(hh_corr, ROOT.TH1D("corr" * l , histtitle * " -- Correlation Layer "  * ll_residuals, 2 * e, -e, e))
        push!(hh_seed_vs_corr, ROOT.TH2D("seed_vs_corr" * l, histtitle * ll_seed_vs_corr, 2 * e, -e, e, 2 * e, -e, e))
        push!(hh_res_vs_bend, ROOT.TH2D("res_vs_bend" * l, histtitle * ll_res_vs_bend, 2 * e, -e, e, 30, -15, 15))
        push!(hh_seed_unbiased, ROOT.TH1D("seed_unbiased" * l , histtitle * " -- Seed Layer "  * ll_pulls, 2 * e, -e, e))
        push!(hh_corr_unbiased, ROOT.TH1D("corr_unbiased" * l , histtitle * " -- Correlation Layer "  * ll_pulls, 2 * e, -e, e))
    end

    h_median = ROOT.TH1D("median", title * ll_residuals, 2 * e, -e, e)
    h_median_unbiased = ROOT.TH1D("median", title * ll_pulls, 100, -3, 3)

    h_mx = ROOT.TH1D("mx", title * "; mx; Events", 100, -0.003, 0.003)
    h_my = ROOT.TH1D("my", title * "; my; Events", 100, -0.003, 0.003)
    h_x0 = ROOT.TH1D("x0", title * "; x0; Events", 100, -5, 5)
    h_y0 = ROOT.TH1D("y0", title * "; x0; Events", 100, -5, 5)

    if ismontecarlo
        h_mx_res = ROOT.TH1D("mx_res", title * "; MC_mx - fit_mx; Events", 100, -0.01, 0.01)
        h_my_res = ROOT.TH1D("my_res", title * "; MC_my - fit_my; Events", 100, -0.01, 0.01)
        h_x0_res = ROOT.TH1D("x0_res", title * "; MC_x0 - fit_x0; Events", 100, -0.01, 0.01)
        h_y0_res = ROOT.TH1D("y0_res", title * "; MC_x0 - fit_x0; Events", 100, -0.01, 0.01)

        h_mx_pulls = ROOT.TH1D("mx_pulls", title * "; (MC_mx - fit_mx)/sigma_mx; Events", 100, -3, 3)
        h_my_pulls = ROOT.TH1D("my_pulls", title * "; (MC_my - fit_my)/sigma_my; Events", 100, -3, 3)
        h_x0_pulls = ROOT.TH1D("x0_pulls", title * "; (MC_x0 - fit_x0)/sigma_x0; Events", 100, -3, 3)
        h_y0_pulls = ROOT.TH1D("y0_pulls", title * "; (MC_x0 - fit_x0)/sigma_x0; Events", 100, -3, 3)
    end

    
end
