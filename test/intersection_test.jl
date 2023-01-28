using Rotations

@testset "intersection" begin
    m = MUonEModule(0.,0.,10.,0.,0.,0.; id=0, name="unknown", spacing=1.8)
    t = Track(0.,0.,0.,0.)
    z = MUonEPede.intersection(m, t)
    mm = modules[2]
    θx, _, _ = Rotations.params(mm.R) 
    zz = MUonEPede.intersection(mm, t)

    @test z ≈ 9.1
    @test zz ≈ mm.r0.z + 0.5 * mm.spacing/cos(θx)
end

@testset "interpolation" begin
    t = Track(-1.1,-1.2,0.1,0.2)
    z = MUonEPede.intersection.(modules, [t])
    s = MUonEPede.local_to_stub.(MUonEPede.global_to_local.(t.(z), modules), modules)
    ss = StubSet(s...)
    tt = MUonEPede.interpolate(ss, modules)
    
    @test t.t0 ≈ tt.t0 rtol=1e-3
    @test t.et ≈ tt.et rtol=1e-3

    s = [Stub(100f0*i, 0.25f0, 0, i) for i in 0:5]
    ss = StubSet(s...)

    tt = MUonEPede.interpolate(ss, modules)

    #@test tt.et.x ≈ 0
    #@test tt.et.y ≈ 0

    z0 = MUonEPede.intersection(modules[1], tt)
    z1 = MUonEPede.intersection(modules[2], tt)
    z4 = MUonEPede.intersection(modules[5], tt)
    z5 = MUonEPede.intersection(modules[6], tt)

    @test MUonEPede.global_to_local(tt(z0), modules[1]).x ≈ MUonEPede.strip_to_local(s[1], modules[1]).x
    @test MUonEPede.global_to_local(tt(z1), modules[2]).x ≈ MUonEPede.strip_to_local(s[2], modules[2]).x
    @test MUonEPede.global_to_local(tt(z4), modules[5]).x ≈ MUonEPede.strip_to_local(s[5], modules[5]).x
    @test MUonEPede.global_to_local(tt(z5), modules[6]).x ≈ MUonEPede.strip_to_local(s[6], modules[6]).x rtol=1.e-3
end
