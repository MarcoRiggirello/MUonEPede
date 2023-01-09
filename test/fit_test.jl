@testset "fit" begin
    t = Track(-1.1,-1.2,0.1,0.2)
    z = intersection.(modules, [t])
    s = MUonEPede.local_to_stub.(global_to_local.(t.(z), modules), modules)
    ss = MUonEPede.StubSet(s...)
    tt, _ = MUonEPede.trackfit(ss, modules)

    @test t.t0 ≈ tt.t0  rtol=1e-3
    @test t.et ≈ tt.et  rtol=1e-3
end