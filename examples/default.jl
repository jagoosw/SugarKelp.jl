using Plots, Dates, Measures

include("../src/Kelp.jl")

t_i = 8.5 * 30
nd = 5 * 365
lat = 45
u = 0.15

#generate the time series as done in paper
u_arr, temp, irr, ex_n = Kelp.defaults(t_i, t_i + nd, u)

a_0 = 0.56;n_0 = 0.018;c_0 = 0.25

#solve
solution, results = Kelp.solvekelp(t_i, nd, u_arr, temp, irr, ex_n, lat, a_0, n_0, c_0)

#plots
pyplot()

display(display(plot(layout=grid(2, 3))))

temp_disp=[];irr_disp=[];n_disp=[];

for t in solution.t
    push!(temp_disp,temp(t))
    push!(irr_disp,irr(t))
    push!(n_disp,ex_n(t))
end

t_disp=solution.t
t_lim=(t_i+nd-365,t_i+nd)

plot!(t_disp,temp_disp,sp=1,ylabel="Temperature/degC",xlim=t_lim)
plot!(t_disp,irr_disp,sp=2,ylabel="Irradiance/micro mol photons / m^2 / s",xlim=t_lim)
plot!(t_disp,n_disp,sp=3,ylabel="Nitrate/micro mol / L",xlim=t_lim)
plot!(t_disp,results.area,sp=4,xlabel="Month", ylabel="Frond Area/dm^2",xlim=t_lim,ylim=(35,45))
plot!(t_disp,results.nitrogen,sp=5, xlabel="Month", ylabel="Nitrogen reserve/gN/g sw",xlim=t_lim)
plot!(t_disp,results.carbon,sp=6, xlabel="Month", ylabel="Carbon reserve/gC/g sw",xlim=t_lim)

xticks=[]
for day in t_i+nd-365:30:nd+t_i
    push!(xticks,Dates.format((Date(2015,8,15)+Dates.Day(day)),"U")[1])
end

plot!(xticks=([t_i+nd-365:30:nd+t_i;],xticks),legend=false,margin=-3mm)

display(display(plot!()))