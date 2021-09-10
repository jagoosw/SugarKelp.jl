using Kelp, Test, Interpolations, DataFrames

@info "These tests check if the code has been broken but does not thoroughly check if it produces correct results (for now at least)"

const arr_lat = [1:3;]
const arr_lon = [1:3;]
const arr_dep = [1:5;]
const arr_t = [1:200;]

const temp = permutedims(repeat([9:(20 - 9) / 199:20;], 1, 5, 4, 3), (4, 3, 2, 1))
const no3 = repeat([5.0], 3, 4, 5, 200)
par = repeat([2.0], 3, 4, 92)
const par_t = [1:2:184;]
const par_fill = -1.0
par[2,3,4] = -1.0
kd = repeat([2.0], 3, 4, 110)
const kd_t = [1:2:220;]
const kd_fill = -1.0
kd[2,3,4] = -1
const beta = repeat([0.1], 3, 4, 5, 200)
const u = repeat([0.15], 3, 4, 5, 200)

t_i = 1.0
nd = 184
a_0 = 0.1;n_0 = 0.022;c_0 = 0.3;

# haveing to use 2013 parameters as cant re define constants but need to be able to run with resp_model 2
@testset "Grid with kd" begin
    all_results = Kelp.solvegrid(t_i, nd, a_0, n_0, c_0, arr_lon, arr_lat, arr_dep, arr_t, no3, temp, u, (par, par_t, par_fill), (kd, kd_t, kd_fill), nothing, "../src/parameters/2013.jl")
    @test size(all_results) == (4, length(arr_lon), length(arr_lat), length(arr_dep), nd + 1)
    @test all_results[1,2,3,4,123] ≈ 0.323811260597751
    @test all_results[2,2,3,4,123] ≈ 0.02159999999999999
    @test all_results[3,2,3,4,123] ≈ 0.010000000000000002
    @test all_results[4,2,3,4,123] ≈ 0.0036635181806166307
end

@testset "Grid with beta" begin
    all_results = Kelp.solvegrid(t_i, nd, a_0, n_0, c_0, arr_lon, arr_lat, arr_dep, arr_t, no3, temp, u, (par, par_t, par_fill), (nothing, nothing, nothing), beta, "../src/parameters/2013.jl")
    @test size(all_results) == (4, length(arr_lon), length(arr_lat), length(arr_dep), nd + 1)
    @test all_results[1,2,3,4,123] ≈ 0.33853815235193496
    @test all_results[2,2,3,4,123] ≈ 0.02159999999999999
    @test all_results[3,2,3,4,123] ≈ 0.010000000000000002
    @test all_results[4,2,3,4,123] ≈ 0.003952014050125887
end

@testset "Grid with beta, specifying resp model to be 1" begin
    all_results = Kelp.solvegrid(t_i, nd, a_0, n_0, c_0, arr_lon, arr_lat, arr_dep, arr_t, no3, temp, u, (par, par_t, par_fill), (nothing, nothing, nothing), beta, "../src/parameters/2013.jl", 1) 
    @test size(all_results) == (4, length(arr_lon), length(arr_lat), length(arr_dep), nd + 1)
    @test all_results[1,2,3,4,123] ≈ 0.33853815235193496
    @test all_results[2,2,3,4,123] ≈ 0.02159999999999999
    @test all_results[3,2,3,4,123] ≈ 0.010000000000000002
    @test all_results[4,2,3,4,123] ≈ 0.003952014050125887
end

@testset "Grid with beta, resp model 2" begin
    all_results = Kelp.solvegrid(t_i, nd, a_0, n_0, c_0, arr_lon, arr_lat, arr_dep, arr_t, no3, temp, u, (par, par_t, par_fill), (nothing, nothing, nothing), beta, "../src/parameters/2013.jl", 2) 
    @test size(all_results) == (4, length(arr_lon), length(arr_lat), length(arr_dep), nd + 1)
    @test all_results[1,2,3,4,123] ≈ 0.6949592729765186
    @test all_results[2,2,3,4,123] ≈ 0.0215977671314431
    @test all_results[3,2,3,4,123] ≈ 0.01000005496349408
    @test all_results[4,2,3,4,123] ≈ 0.01014604344770969
end 

@testset "Grid with beta, resp model 2, specifying dt=1" begin
    all_results = Kelp.solvegrid(t_i, nd, a_0, n_0, c_0, arr_lon, arr_lat, arr_dep, arr_t, no3, temp, u, (par, par_t, par_fill), (nothing, nothing, nothing), beta, "../src/parameters/2013.jl", 1, 1) 
    @test size(all_results) == (4, length(arr_lon), length(arr_lat), length(arr_dep), nd + 1)
    @test all_results[1,2,3,4,123] ≈ 0.33853815235193496
    @test all_results[2,2,3,4,123] ≈ 0.02159999999999999
    @test all_results[3,2,3,4,123] ≈ 0.010000000000000002
    @test all_results[4,2,3,4,123] ≈ 0.003952014050125887
end
@testset "Grid with beta, resp model 2, specifying dt=2" begin
    # This configuration is unstable so area goes below zero and it stops, that would be a useful behaviour to test for but its not done here
    all_results = Kelp.solvegrid(t_i, nd, a_0, n_0, c_0, arr_lon, arr_lat, arr_dep, arr_t, no3, temp, u, (par, par_t, par_fill), (nothing, nothing, nothing), beta, "../src/parameters/2013.jl", 1, 2) 
    @test size(all_results) == (4, length(arr_lon), length(arr_lat), length(arr_dep), round(nd / 2 + 1))
    @test all_results[1,2,3,4,5] ≈ 0.19669013321363188
    @test all_results[2,2,3,4,5] ≈ 0.01614271696738769
    @test all_results[3,2,3,4,5] ≈ 0.02894612990034886
    @test all_results[4,2,3,4,5] ≈ 0.0011375879549344894
end
@testset "Grid with beta, resp model 2, specifying dt=1, specifying to show progress" begin
    all_results = Kelp.solvegrid(t_i, nd, a_0, n_0, c_0, arr_lon, arr_lat, arr_dep, arr_t, no3, temp, u, (par, par_t, par_fill), (nothing, nothing, nothing), beta, "../src/parameters/2013.jl", 1, 1, true) 
    @test size(all_results) == (4, length(arr_lon), length(arr_lat), length(arr_dep), nd + 1)
    @test all_results[1,2,3,4,123] ≈ 0.33853815235193496
    @test all_results[2,2,3,4,123] ≈ 0.02159999999999999
    @test all_results[3,2,3,4,123] ≈ 0.010000000000000002
    @test all_results[4,2,3,4,123] ≈ 0.003952014050125887
end
# This is a test for progress supression but there is no way to actually test it
# all_results = Kelp.solvegrid(t_i, nd, a_0, n_0, c_0, arr_lon, arr_lat, arr_dep, arr_t, no3, temp, u, (par, par_t, par_fill), (nothing, nothing, nothing), beta, "../src/parameters/2013.jl", 1, 1, false)

u_itp = Interpolations.LinearInterpolation(arr_t, u[1,1,1,:])
temp_itp = Interpolations.LinearInterpolation(arr_t, temp[1,1,1,:])
no3_itp = Interpolations.LinearInterpolation(arr_t, no3[1,1,1,:])
par_itp = Interpolations.LinearInterpolation(par_t, par[1,1,:], extrapolation_bc=Flat())

@testset "Single point, dataframe output" begin 
    solution, a_result = Kelp.solvekelp(t_i, nd, u_itp, temp_itp, par_itp, no3_itp, 55, a_0, n_0, c_0, "../src/parameters/2013.jl")
    @test typeof(a_result) == DataFrames.DataFrame
    @test size(a_result) == (185, 5)
    @test a_result.area[17] ≈ 0.36857308566331204
    @test a_result.nitrogen[17] ≈ 0.01579264171780362
    @test a_result.carbon[17] ≈ 0.027214793958727508
end

@testset "Depth interpolator utility" begin
    tnew = Kelp.interp_deps(temp, arr_dep, [1:0.5:5;], -1)
    @test size(tnew) == (3, 4, 9, 200)
    @test tnew[2,3,4,5] ≈ 9.22110552763819
end

@testset "Indexer" begin @test Kelp.get_ind(2, arr_lat, 0) == 2; @test Kelp.get_ind(2.34, arr_lat, 1) == 2 end