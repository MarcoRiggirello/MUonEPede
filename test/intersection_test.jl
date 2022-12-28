using Rotations

@testset "intersection" begin
    m = MUonEModule(0.,0.,10.,0.,0.,0.; id=0, name="unknown", spacing=1.8)
    t = Track(0.,0.,0.,0.)
    z = intersection(m, t)
    mm = modules[2]
    θx, _, _ = Rotations.params(mm.R) 
    zz = intersection(mm, t)

    @test z ≈ 9.1
    @test zz ≈ mm.r0.z + 0.5 * mm.spacing/cos(θx)
end

@testset "interpolation" begin
    t = Track(1.1,2.2,0.1,0.2)
    h0 = t(10)
    h1 = t(20)
    h4 = t(50)
    h5 = t(60)
    tt = interpolate(h0, h1, h4, h5)
    
    @test t.t0 ≈ tt.t0 
    @test t.et ≈ tt.et 
    @test tt.et.x ≈ 0.1
    @test tt.et.y ≈ 0.2

    l0 = strip_to_local(100f0, 0.25f0, modules[1])
    l1 = strip_to_local(200f0, 0.75f0, modules[2])
    l4 = strip_to_local(100f0, 0.25f0, modules[5])
    l5 = strip_to_local(200f0, 0.75f0, modules[6])

    tt = interpolate(
        local_to_global(l0, modules[1]),
        local_to_global(l1, modules[2]),
        local_to_global(l4, modules[5]),
        local_to_global(l5, modules[6])
    )

    @test tt.et.x ≈ 0
    @test tt.et.y ≈ 0

    z0 = intersection(modules[1], tt)
    z1 = intersection(modules[2], tt)
    z4 = intersection(modules[5], tt)
    z5 = intersection(modules[6], tt)

    @test global_to_local(tt(z0), modules[1]).x ≈ l0.x
    @test global_to_local(tt(z1), modules[2]).x ≈ l1.x
    @test global_to_local(tt(z4), modules[5]).x ≈ l4.x
    @test global_to_local(tt(z5), modules[6]).x ≈ l5.x
end
