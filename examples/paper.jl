using Plots, Dates, Measures, Interpolations, RollingFunctions, NumericalIntegration, CSV, DataFrames

include("../src/Kelp.jl")
include("../src/parameters.jl")

offset = 4*31+28+30*2.0
# collected from fig2 of origional paper with https://www.graphreader.com/
temp_file=CSV.read("examples/bs2012/temp.csv",DataFrame);sort!(temp_file);temp_t=temp_file.day;temp=temp_file.temp
no3_file=CSV.read("examples/bs2012/no3.csv",DataFrame);sort!(no3_file);no3_t=no3_file.day;no3=no3_file.no3
irr_file=CSV.read("examples/bs2012/irr.csv",DataFrame);sort!(irr_file);irr_t=irr_file.day;irr=irr_file.irr
temp_t .+= offset ;no3_t .+= offset ;irr_t .+= offset

#irr .*= exp(-0.07*13) # Background attenuation is set to 0.07m^-1(?)

t_i = offset+12
nd = 365+16
lat = 60.257
u = 0.15

u_arr = Interpolations.LinearInterpolation([t_i:t_i + nd;], fill(u, Int(nd + 1)))

temp_arr = Interpolations.LinearInterpolation(temp_t, temp, extrapolation_bc=Flat())
no3_arr = Interpolations.LinearInterpolation(no3_t, no3, extrapolation_bc=Flat())
irr_arr = Interpolations.LinearInterpolation(irr_t, irr, extrapolation_bc=Flat())

a_0 = 30;n_0 = 0.01;c_0 = 0.6

solution, results = Kelp.solvekelp(t_i, nd, u_arr, temp_arr, irr_arr, no3_arr, lat, a_0, n_0, c_0)

pyplot()
plot(layout=grid(2, 3))

plot!(temp_t,temp,sp=1,ylabel="Temperature/degC",legend=true)
plot!(irr_t,irr,sp=2,ylabel="Irradiance/micro mol photons / m^2 / s",legend=true)
plot!(no3_t,no3,sp=3,ylabel="Nitrate/micro mol / L",legend=true)

# sjotun 1993
c_t = [10.016,100.156,130.973,164.872,192.607,224.965,266.568,291.222,343.611,409.097] .+ offset 
c = [0.312,0.349,0.292,0.195,0.224,0.266,0.281,0.307,0.363,0.373]
n_t = [10.255,96.239,135.681,160.924,195.633,222.454,269.785,288.717,345.514,406.255] .+ offset 
n = [0.007,0.014,0.021,0.024,0.025,0.022,0.023,0.017,0.012,0.013]
area_file=CSV.read("examples/bs2012/area.csv",DataFrame);sort!(area_file);area_t=area_file.day .+offset;area=area_file.area
c2_file=CSV.read("examples/bs2012/carbon.csv",DataFrame);sort!(c2_file);c2_t=c2_file.day .+offset;c2=c2_file.carbon
n2_t = [14.8696,32.8696,57.913,75.1304,90.7826,104.8696,125.2174,145.5652,156.5217,162.7826,170.6087,178.4348,184.6957,194.087,200.3478,214.4348,237.913,248.8696,259.8261,267.6522,270.7826,273.913,278.6087,288,297.3913,308.3478,316.1739,327.1304,336.5217,352.1739,364.6957,378.7826,396.7826] .+ offset 
n2 = [0.0088,0.0095,0.0109,0.0124,0.0146,0.017,0.0209,0.0249,0.027,0.0284,0.0291,0.0285,0.0281,0.0265,0.025,0.0229,0.0189,0.0175,0.0165,0.0154,0.0145,0.0136,0.0127,0.0112,0.0102,0.0095,0.0091,0.0086,0.0083,0.0083,0.0084,0.0086,0.0088]
# structural to dry weight conversion (paper plots g/g dry weight where as g/g structural weight is used in calculations), I think they are also plotting the total carbon not just the reserve
n_factor = (results.nitrogen .- N_min) .* K_N
c_factor = (results.carbon .- C_min) .* K_C
w_factor = 1 .+ n_factor .+ c_factor .+ C_min .+ N_min

plot!(results.time,results.area,sp=4,xlabel="Month", ylabel="Frond Area/dm^2",ylim=(29, 46),label="Model")
plot!(area_t,area,label="Broch and Slagstad, 2012 (model)",sp=4)

plot!(results.time,(results.nitrogen .+ N_struct) ./ w_factor,sp=5, xlabel="Month", ylabel="Nitrogen reserve/gN/g dw",label="Model")
plot!(n_t,n,sp=5,seriestype=:scatter, label="Sjotun, 1993")
plot!(n2_t,n2,label="Broch and Slagstad, 2012 (model)",sp=5)

plot!(results.time,(results.carbon .+ C_struct) ./ w_factor,sp=6, xlabel="Month", ylabel="Carbon reserve/gC/g dw",label="Model")
plot!(c_t,c,sp=6,seriestype=:scatter, label="Sjotun, 1993",legend=:bottomleft)
plot!(c2_t,c2,label="Broch and Slagstad, 2012 (model)",sp=6)

t_ticks = []
val_ticks = []
for day in t_i:t_i + nd
    date = Date(1981, 1, 1) + Dates.Day(day)
    if Dates.format(date, "d") == "1"
        push!(t_ticks, day)
        push!(val_ticks, Dates.format(date, "U")[1])
    end
end

plot!(xticks=(t_ticks, val_ticks),margin=-3mm)
plot!(top_margin=0mm,legend=false)
display(display(plot!()))
