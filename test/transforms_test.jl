@testset "transforms" begin
    mu = MUonEModule(0,0,12,0,0.233,pi; name="dummyU", id=0, spacing=4.0)
    my = MUonEModule(0,0,0,-0.233,pi,pi/2; name="dummyY", id=1, spacing=1.8)
    lu = strip_to_local(Stub(500, 0.25, 0, 0), mu)
    ly = strip_to_local(Stub(100, 0.75, 0, 1), my)
    
    @test lu.z ≈ -2.0
    @test ly.z ≈ -0.9
    @test lu ≈ global_to_local(local_to_global(lu, mu), mu)
    @test ly ≈ global_to_local(local_to_global(ly, my), my)
    @test local_to_global(ly, my).y ≈ ly[1] * cos(0.233) - ly[3] * sin(0.233) 
end
