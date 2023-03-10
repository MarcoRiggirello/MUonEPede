@testset "transforms" begin
    mu = MUonEModule(0,0,12,0,0.233,pi; name="dummyU", id=0, spacing=4.0)
    my = MUonEModule(0,0,0,pi,0.233,-pi/2; name="dummyY", id=1, spacing=1.8)
    lu = MUonEPede.stub_to_local(Stub(500, 0.25, 0, 0), mu)[1]
    ly = MUonEPede.stub_to_local(Stub(100, 0.75, 0, 1), my)[1]
    lc = MUonEPede.stub_to_local(Stub(100, 0.75, 0, 1), my)[2]
    
    @test lu.z ≈ -2.0
    @test ly.z ≈ -0.9
    @test lc.z ≈ +0.9
    @test lu ≈ MUonEPede.global_to_local(MUonEPede.local_to_global(lu, mu), mu)
    @test ly ≈ MUonEPede.global_to_local(MUonEPede.local_to_global(ly, my), my)
    #@test MUonEPede.local_to_global(ly, my).y ≈ ly[1] * cos(0.233) - ly[3] * sin(0.233) 

    lc0 = @SVector [0.0, -2.5, 0.0] # CIC 0 
    lc1 = @SVector [0.0, 2.5, 0.0] # CIC 1
    lss = @SVector [-1.0, 0.0, 0.0] # SEH side
    lso = @SVector [1.0, 0.0, 0.0] # opposite to SEH
    
    lt = @SVector [0.0, 0.0, 0.09] # top sensor
    lb = @SVector [0.0, 0.0, -0.09] # bottom sensor

    # X module
    @test MUonEPede.local_to_global(lb, modules[1]).z < MUonEPede.local_to_global(lt, modules[1]).z 
    @test MUonEPede.local_to_global(lss, modules[1]).x > 0 
    @test MUonEPede.local_to_global(lso, modules[1]).x < 0
    @test MUonEPede.local_to_global(lc0, modules[1]).y > 0 
    @test MUonEPede.local_to_global(lc1, modules[1]).y < 0
    # Y module
    @test MUonEPede.local_to_global(lb, modules[2]).z > MUonEPede.local_to_global(lt, modules[2]).z 
    @test MUonEPede.local_to_global(lss, modules[2]).y < 0 
    @test MUonEPede.local_to_global(lso, modules[2]).y > 0
    @test MUonEPede.local_to_global(lc0, modules[2]).x < 0 
    @test MUonEPede.local_to_global(lc1, modules[2]).x > 0
    # U module
    @test MUonEPede.local_to_global(lb, modules[3]).z > MUonEPede.local_to_global(lt, modules[3]).z 
    @test MUonEPede.local_to_global(lss, modules[3]).x > 0 
    @test MUonEPede.local_to_global(lss, modules[3]).y < 0 
    @test MUonEPede.local_to_global(lso, modules[3]).x < 0 
    @test MUonEPede.local_to_global(lso, modules[3]).y > 0 
    @test MUonEPede.local_to_global(lc0, modules[3]).x < 0
    @test MUonEPede.local_to_global(lc0, modules[3]).y < 0
    @test MUonEPede.local_to_global(lc1, modules[3]).x > 0
    @test MUonEPede.local_to_global(lc1, modules[3]).y > 0
    @test MUonEPede.local_to_global(lss, modules[3]).x ≈ -MUonEPede.local_to_global(lss, modules[3]).y rtol=1.e-4 
    # V module
    @test MUonEPede.local_to_global(lb, modules[4]).z < MUonEPede.local_to_global(lt, modules[4]).z 
    @test MUonEPede.local_to_global(lss, modules[4]).x < 0 
    @test MUonEPede.local_to_global(lss, modules[4]).y < 0 
    @test MUonEPede.local_to_global(lso, modules[4]).x > 0 
    @test MUonEPede.local_to_global(lso, modules[4]).y > 0 
    @test MUonEPede.local_to_global(lc0, modules[4]).x > 0
    @test MUonEPede.local_to_global(lc0, modules[4]).y < 0
    @test MUonEPede.local_to_global(lc1, modules[4]).x < 0
    @test MUonEPede.local_to_global(lc1, modules[4]).y > 0
    @test MUonEPede.local_to_global(lss, modules[4]).x ≈ MUonEPede.local_to_global(lss, modules[4]).y rtol=1.e-4
    # X module
    @test MUonEPede.local_to_global(lb, modules[5]).z < MUonEPede.local_to_global(lt, modules[5]).z 
    @test MUonEPede.local_to_global(lss, modules[5]).x > 0 
    @test MUonEPede.local_to_global(lso, modules[5]).x < 0
    @test MUonEPede.local_to_global(lc0, modules[5]).y > 0 
    @test MUonEPede.local_to_global(lc1, modules[5]).y < 0
    # Y module
    @test MUonEPede.local_to_global(lb, modules[6]).z > MUonEPede.local_to_global(lt, modules[6]).z 
    @test MUonEPede.local_to_global(lss, modules[6]).y < 0 
    @test MUonEPede.local_to_global(lso, modules[6]).y > 0
    @test MUonEPede.local_to_global(lc0, modules[6]).x < 0 
    @test MUonEPede.local_to_global(lc1, modules[6]).x > 0
end
