using Plots

include("../src/Kelp.jl")

t_i = 4.5 * 30
nd = 5 * 365
lat = 45
u = 0.15

u, temp, irr, ex_n = Kelp.defaults(t_i, t_i + nd, u)

a_0 = 0.56;n_0 = 0.018;c_0 = 0.25

solution, results = Kelp.solvekelp(t_i, nd, u, temp, irr, ex_n, lat, a_0, n_0, c_0)


pyplot()

display(display(plot(layout=grid(2, 2))))
plot!(solution.t,results.area,sp=1,title="Frond Area",xlim=(solution.t[length(solution.t)] - 365, solution.t[length(solution.t)]), xlabel="Time/days", ylabel="Frond Area/dm^2")
plot!(solution.t,results.nitrogen,sp=2,title="Nitrogen Content",xlim=(solution.t[length(solution.t)] - 365, solution.t[length(solution.t)]), xlabel="Time/days", ylabel="Nitrogen reserve/gN/g sw")
plot!(solution.t,results.carbon,sp=3,title="Carbon Content",xlim=(solution.t[length(solution.t)] - 365, solution.t[length(solution.t)]), xlabel="Time/days", ylabel="Carbon reserve/gC/g sw")
display(display(plot!()))