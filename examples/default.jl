using Plots, Dates, Measures

include("../src/Kelp.jl")

t_i=8.0*30
nd=365
lat=60.25
u=0.15

#generate the time series as done in paper
u_arr, temp, irr, ex_n = Kelp.defaults(t_i, t_i + nd, u)

a_0 = 30;n_0 = 0.009;c_0 = 0.3

#solve
solution, results = Kelp.solvekelp(t_i, nd, u_arr, temp, irr, ex_n, lat, a_0, n_0, c_0)

#plots
pyplot()

temp_disp=[];irr_disp=[];n_disp=[];

for t in solution.t
    push!(temp_disp,temp(t))
    push!(irr_disp,irr(t))
    push!(n_disp,ex_n(t))
end

t_disp=solution.t
t_lim=(t_i,t_i+365)

display(display(plot(layout=grid(2, 3))))

plot!(t_disp,temp_disp,sp=1,ylabel="Temperature/degC",xlim=t_lim,label="Analytical inputs")
plot!(t_disp,irr_disp,sp=2,ylabel="Irradiance/micro mol photons / m^2 / s",xlim=t_lim,label="Analytical inputs")
plot!(t_disp,n_disp,sp=3,ylabel="Nitrate/micro mol / L",xlim=t_lim,label="Analytical inputs")
plot!(t_disp,results.area,sp=4,xlabel="Month", ylabel="Frond Area/dm^2",xlim=t_lim,ylim=(30,45),label="Analytical inputs")
plot!(t_disp,results.nitrogen,sp=5, xlabel="Month", ylabel="Nitrogen reserve/gN/g sw",xlim=t_lim,label="Analytical inputs")
plot!(t_disp,results.carbon,sp=6, xlabel="Month", ylabel="Carbon reserve/gC/g sw",xlim=t_lim,label="Analytical inputs")

t_ticks = []
val_ticks = []
for day in t_i:t_i + nd
    date=Date(1981, 1, 1) + Dates.Day(day)
    if Dates.format(date,"d")=="1"
        push!(t_ticks, day)
        push!(val_ticks, Dates.format(date, "U")[1])
    end
end

plot!(xticks=(t_ticks, val_ticks),margin=-3mm)#legend=false)

display(display(plot!()))