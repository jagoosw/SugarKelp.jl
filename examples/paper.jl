using Plots, Dates, Measures, Interpolations, RollingFunctions, NumericalIntegration, CSV, DataFrames

include("../src/Kelp.jl")
include("../src/parameters.jl")

offset = 4*31+28+30*2.0
# collected from fig2 of origional paper with https://www.graphreader.com/
temp_file=CSV.read("examples/b&s_temp.csv",DataFrame);sort!(temp_file);temp_t=temp_file.day;temp=temp_file.temp
no3_file=CSV.read("examples/b&s_no3.csv",DataFrame);sort!(no3_file);no3_t=no3_file.day;no3=no3_file.no3
irr_file=CSV.read("examples/b&s_irr.csv",DataFrame);sort!(irr_file);irr_t=irr_file.day;irr=irr_file.irr

# no3 ./= 1
irr .*= 10^6 / (24 * 60 * 60)#
temp_t .+= offset ;no3_t .+= offset ;irr_t .+= offset

irr .*= exp(-0.07*13) # Background attenuation is set to 0.07m^-1(?)

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

t_disp = solution.t

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
area_file=CSV.read("examples/b&s_area.csv",DataFrame);sort!(area_file);area_t=area_file.day .+offset;area=area_file.area
c2_file=CSV.read("examples/b&s_carbon.csv",DataFrame);sort!(c2_file);c2_t=c2_file.day .+offset;c2=c2_file.carbon
n2_t = [14.8696,32.8696,57.913,75.1304,90.7826,104.8696,125.2174,145.5652,156.5217,162.7826,170.6087,178.4348,184.6957,194.087,200.3478,214.4348,237.913,248.8696,259.8261,267.6522,270.7826,273.913,278.6087,288,297.3913,308.3478,316.1739,327.1304,336.5217,352.1739,364.6957,378.7826,396.7826] .+ offset 
n2 = [0.0088,0.0095,0.0109,0.0124,0.0146,0.017,0.0209,0.0249,0.027,0.0284,0.0291,0.0285,0.0281,0.0265,0.025,0.0229,0.0189,0.0175,0.0165,0.0154,0.0145,0.0136,0.0127,0.0112,0.0102,0.0095,0.0091,0.0086,0.0083,0.0083,0.0084,0.0086,0.0088]
# structural to dry weight conversion (paper plots g/g dry weight where as g/g structural weight is used in calculations), I think they are also plotting the total carbon not just the reserve
n_factor = (results.nitrogen .- N_min) .* K_N
c_factor = (results.carbon .- C_min) .* K_C
w_factor = 1 .+ n_factor .+ c_factor .+ C_min .+ N_min

plot!(t_disp,results.area,sp=4,xlabel="Month", ylabel="Frond Area/dm^2",ylim=(30, 45),label="Model")
plot!(area_t,area,label="Paper",sp=4)

plot!(t_disp,(results.nitrogen .+ N_struct) ./ w_factor,sp=5, xlabel="Month", ylabel="Nitrogen reserve/gN/g dw",label="Model")
plot!(n_t,n,sp=5,seriestype=:scatter, label="Sjotun, 1993")
plot!(n2_t,n2,label="Paper",sp=5)

plot!(t_disp,(results.carbon .+ C_struct) ./ w_factor,sp=6, xlabel="Month", ylabel="Carbon reserve/gC/g dw",label="Model")
plot!(c_t,c,sp=6,seriestype=:scatter, label="Sjotun, 1993",legend=:bottomleft)
plot!(c2_t,c2,label="Paper",sp=6)

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
plot!(top_margin=0mm)
display(display(plot!()))

c_fixed = integrate(solution.t[2:end], diff(results.carbon) .* (K_A .* results.area[2:end]))

mu = []

delts = []
for j = 1:365
    theta = 0.2163108 + 2 * atan(0.9671396 * tan(0.00860 * (j - 186))) # revolution angle from day of the year
    dec = asin(0.39795 * cos(theta)) # sun declination angle 
    p = 0.8333 # sunrise/sunset is when the top of the sun is apparently even with horizon
    push!(
            delts,
            24 -
            (24 / pi) * acos(
                (sin(p * pi / 180) + sin(lat * pi / 180) * sin(dec)) /
                (cos(lat * pi / 180) * cos(dec)),
            ),
        )
end
DeltaL = diff(delts)
push!(DeltaL, DeltaL[364])
NormDeltaL = DeltaL / findmax(DeltaL)[1]

    
for (ind, t) in enumerate(solution.t)
    d = trunc(Int, mod(floor(t), 365) + 1)
    lambda =  NormDeltaL[d] # This seems wrong to be interpolated because "change in day length" is discrete so going to stick choosing day numbers
    # On the other hand the other parameters could be time series with higher resolution than once per day (but again is that valid with this model because they just one per day)

    f_area = m_1 * exp(-(results.area[ind] / A_0)^2) + m_2 # effect of size on growth rate
    
    temp = temp_arr(t)
    
    if temp < 10 # effect of temperature on growth rate
        f_temp = 0.08 * temp + 0.2
    elseif temp < 15
        f_temp = 1
    elseif temp < 19
        f_temp = 19 / 4 - temp / 4
    else
        f_temp = 0
    end
    
    f_photo = a_1 * (1 + sign(lambda) * sqrt(abs(lambda))) + a_2
        
    push!(mu, f_area * f_temp * f_photo * min(1 - N_min / results.nitrogen[ind], 1 - C_min / results.carbon[ind]))
end

gross_a = integrate(solution.t, mu .* results.area)

println("Total carbon fixed: $c_fixed /g")
println("Gross frond area: $gross_a /dm^2")