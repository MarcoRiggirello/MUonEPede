@testset "transforms" begin
    mu = MUonEModule(0,0,12,0,pi,pi/4; name="dummyX", id=0, type='U')
    my = MUonEModule(0,0,30,0.233,0,0; name="dummyY", id=1, type='Y')
    lu = MUonEPede.strip_to_local(1000, 0.25, mu)
    ly = MUonEPede.strip_to_local(100, 0.75, my)
    
    @test lu.z ≈ -2.0
    @test ly.z ≈ -0.9
    @test lu ≈ global_to_local(local_to_global(lu, mu), mu)
    @test ly ≈ global_to_local(local_to_global(ly, my), my)
    @test local_to_global(ly, my).y ≈ ly[1] * cos(0.233)
end
