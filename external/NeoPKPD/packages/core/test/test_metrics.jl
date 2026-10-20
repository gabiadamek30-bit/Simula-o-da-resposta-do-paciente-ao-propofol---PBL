# PK/PD metrics tests

@testset "Exposure metrics: trapezoid AUC and Cmax" begin
    t = [0.0, 1.0, 2.0]
    c = [0.0, 1.0, 1.0]
    @test cmax(t, c) == 1.0
    # AUC = 0.5*(0+1)*1 + 0.5*(1+1)*1 = 0.5 + 1.0 = 1.5
    @test isapprox(auc_trapezoid(t, c), 1.5; rtol = 0.0, atol = 1e-12)
end

@testset "Response metrics: emin, time_below, auc_above_baseline" begin
    t = [0.0, 1.0, 2.0]
    y = [10.0, 8.0, 9.0]
    @test emin(t, y) == 8.0
    # left-constant rule: interval [0,1] uses y[1]=10 >= 9.5 (not counted),
    # interval [1,2] uses y[2]=8 < 9.5 (counted), so total = 1.0
    @test isapprox(time_below(t, y, 9.5), 1.0; atol = 1e-12)
    # baseline 10: suppression curve is [0,2,1]; AUC = 0.5*(0+2)*1 + 0.5*(2+1)*1 = 1 + 1.5 = 2.5
    @test isapprox(auc_above_baseline(t, y, 10.0), 2.5; atol = 1e-12)
end
