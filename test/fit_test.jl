using Rotations: Random
@testset "fit" begin
    t = Track(-1.1f0,-1.2f0,0.001f0,0.002f0)
    z = MUonEPede.intersection.(modules, [t])
    z_s = getindex.(z, 1)
    z_c = getindex.(z, 2)
    ρ_s = MUonEPede.global_to_local.(t.(z_s), modules)
    ρ_c = MUonEPede.global_to_local.(t.(z_c), modules)
    s = MUonEPede.local_to_stub.(ρ_s, ρ_c, modules, 0.0f0)
    ss = StubSet(s...)
    tt = MUonEPede.trackfit(ss, modules, 1.0f0)

    @test t.t0 ≈ tt.t0  rtol=1#e-3
    @test t.et ≈ tt.et  rtol=1#e-3
    @show tt.t0
    @show tt.et

    p_0 = [40.f0,50.f0,10.f0,-20.f0]

    fg!(F, G, x) = MUonEPede.chisquare_gradient!(F, G, ss, modules, Track(x...), 1.0f0)
    results = optimize(Optim.only_fg!(fg!), p_0, BFGS(linesearch = BackTracking()))
    popt = Optim.minimizer(results)
    ttt = Track(popt...)

    #@test ttt.t0 ≈ t.t0  rtol=1e-3
    #@test ttt.et ≈ t.et  rtol=1e-3
end
