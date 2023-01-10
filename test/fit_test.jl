@testset "fit" begin
    t = Track(-1.1,-1.2,0.1,0.2)
    z = MUonEPede.intersection.(modules, [t])
    s = MUonEPede.local_to_stub.(MUonEPede.global_to_local.(t.(z), modules), modules)
    ss = StubSet(s...)
    tt, _ = MUonEPede.trackfit(ss, modules)

    @test t.t0 ≈ tt.t0  rtol=1e-3
    @test t.et ≈ tt.et  rtol=1e-3
end