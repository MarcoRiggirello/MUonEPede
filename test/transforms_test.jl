@testset "transforms" begin
    mu = MUonEModule(0,0,12,0,0.233,pi; name="dummyU", id=0, spacing=4.0)
    my = MUonEModule(0,0,0,-0.233,pi,pi/2; name="dummyY", id=1, spacing=1.8)
    lu = MUonEPede.strip_to_local(Stub(500, 0.25, 0, 0), mu)
    ly = MUonEPede.strip_to_local(Stub(100, 0.75, 0, 1), my)
    
    @test lu.z ≈ -2.0
    @test ly.z ≈ -0.9
    @test lu ≈ MUonEPede.global_to_local(MUonEPede.local_to_global(lu, mu), mu)
    @test ly ≈ MUonEPede.global_to_local(MUonEPede.local_to_global(ly, my), my)
    @test MUonEPede.local_to_global(ly, my).y ≈ ly[1] * cos(0.233) - ly[3] * sin(0.233) 

    l = @SVector [3.0, 2.5, -0.09]

    # X module
    @test MUonEPede.local_to_global(l, modules[1]).x < 0 
    @test MUonEPede.local_to_global(l, modules[1]).y < 0
    # Y module
    @test MUonEPede.local_to_global(l, modules[2]).x > 0 
    @test MUonEPede.local_to_global(l, modules[2]).y > 0
    # U module
    @test MUonEPede.local_to_global(l, modules[3]).x < 0 
    @test MUonEPede.local_to_global(l, modules[3]).y > 0
    # V module
    @test MUonEPede.local_to_global(l, modules[4]).x > 0 
    @test MUonEPede.local_to_global(l, modules[4]).y > 0
end
