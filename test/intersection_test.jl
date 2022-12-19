@testset "interpolation" begin
    t = Track(1.1,2.2,0.1,0.2)
    h0 = t(10)
    h1 = t(20)
    h4 = t(50)
    h5 = t(60)
    tt = interpolate(h0, h1, h4, h5)
    
    @test t.t0 ≈ tt.t0 
    @test t.et ≈ tt.et 
    @test tt.θ ≈ 0.1
    @test tt.ϕ ≈ 0.2

end