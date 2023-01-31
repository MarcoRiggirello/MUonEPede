using Rotations: Random
@testset "fit" begin
    t = Track(-1.1,-1.2,0.1,0.2)
    z = getindex.(MUonEPede.intersection.(modules, [t]),1)
    s = MUonEPede.local_to_stub.(MUonEPede.global_to_local.(t.(z), modules), modules)
    ss = StubSet(s...)
    tt, _ = MUonEPede.trackfit(ss, modules, 1.0)

    @test t.t0 ≈ tt.t0  rtol=1e-3
    @test t.et ≈ tt.et  rtol=1e-3

    p_0 = [40.,50.,10.,-20.]

    fg!(F, G, x) = MUonEPede.chisquare_gradient!(F, G, ss, modules, Track(x...))
    results = optimize(Optim.only_fg!(fg!), p_0, BFGS(linesearch = BackTracking()))
    popt = Optim.minimizer(results)
    ttt = Track(popt...)

    @test ttt.t0 ≈ t.t0  rtol=1e-3
    @test ttt.et ≈ t.et  rtol=1e-3
end
