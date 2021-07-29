module Kelp
using RecursiveArrayTools, DiffEqBase, OrdinaryDiffEq, Roots, Interpolations, DataFrames
# Input parameters
# t: time (days)
# irr: irradiance (micro mol photons / m^2 / s)
# temp: Temperature (deg C)
# ex_n: External (to kelp) nitrate concentration (micro mol / L) 
# u: Current speed (relative to kelp) (m/s)

# Variables
# a: Frond area of individual kelp (dm^2)
# n: Nitrogen reserve relative to dry weight (gN/(g sw))
# c: Carbon reserve relative to dry weight (gC/(g sw))
function equations!(y, params, t)
    a, n, c = y

    if a > 0
        if c <= C_min 
            a_check = a - a * (C_min - c) / C_struct
        else
            a_check = 0
        end

        u_arr, temp_arr, irr_arr, ex_n_arr, NormDeltaL, resp_model, dt = params

        u = u_arr(t)
        temp = temp_arr(t)
        irr = irr_arr(t)
        ex_n = ex_n_arr(t)

        p_max =
                P_1 * exp(T_AP / T_P1 - T_AP / (temp + 273.15)) / (
                    1 +
                    exp(T_APL / (temp + 273.15) - T_APL / T_PL) +
                    exp(T_APH / T_PH - T_APH / (temp + 273.15))
                )

        beta_func(x) = p_max - (alpha * I_sat / log(1 + alpha / x)) * (alpha / (alpha + x)) * (x / (alpha + x))^(x / alpha)
        beta = find_zero(beta_func, (0, 0.1), Bisection())

        p_s = alpha * I_sat / log(1 + alpha / beta)
        p = p_s * (1 - exp(-alpha * irr / p_s)) * exp(-beta * irr / p_s) # gross photosynthesis
        e = 1 - exp(gamma * (C_min - c)) # carbon exudation

        d = trunc(Int, mod(floor(t), 365) + 1) 
        lambda = NormDeltaL[d] 
        
        f_area = m_1 * exp(-(a / A_0)^2) + m_2 # effect of size on growth rate

        if -1.8 <= temp < 10 # effect of temperature on growth rate
            f_temp = 0.08 * temp + 0.2
        elseif 10 <= temp <= 15
            f_temp = 1
        elseif 15 < temp <= 19
            f_temp = 19 / 4 - temp / 4
        elseif temp > 19
            f_temp = 0
        else
            throw("Out of range temperature, $temp")
        end

        f_photo = a_1 * (1 + sign(lambda) * abs(lambda)^.5) + a_2

        mu = f_area * f_temp * f_photo * min(1 - N_min / n, 1 - C_min / c) # gross area specific growth rate
        nu = 1e-6 * exp(epsilon * a) / (1 + 1e-6 * (exp(epsilon * a) - 1)) # front erosion
        j = J_max * (ex_n / (K_X + ex_n)) * ((N_max - n) / (N_max - N_min)) * (1 - exp(-u / U_0p65)) # nitrate uptake rate

        if resp_model == 1
            r = R_1 * exp(T_AR / T_R1 - T_AR / (temp + 273.15)) # temperature dependent respiration
        elseif resp_model == 2
            r = (R_A * (mu / mu_max + j / J_max) + R_B) * exp(T_AR / T_R1 - T_AR / (temp + 273.15))
        end

        da = (mu - nu) * a
        dn = j / K_A - mu * (n + N_struct)
        dc = (p * (1 - e) - r) / K_A - mu * (c + C_struct)
            
        if c + dc < C_min
            da -= (a * (C_min - c) / C_struct) * .5 * dt
            dc = (C_min - c) * .5 * dt
        end
    else
        da, dn, dc, j = 0, 0, 0,0 
        @warn "Area reached 0"
    end
    return (vcat(da, dn, dc, j*a))
end

function defaults(t_i, t_e, u)
    t_arr, irr_arr, ex_n_arr = [], [], []
    for t = t_i:t_e
        d = trunc(Int, mod(floor(t), 365) + 1)
        push!(t_arr, 6 * cos((d - 250) * 2 * pi / 365) + 8)
        push!(irr_arr, 40 * (sin((d + 15) * pi / 365)^10) + 1)
        push!(ex_n_arr, (7 * ((cos(d * 2 * pi / 365) + 1) / 2).^3 + 0.1) / 1000)
    end

    t_arr = Interpolations.LinearInterpolation([t_i:t_e;], t_arr)
    irr_arr = Interpolations.LinearInterpolation([t_i:t_e;], irr_arr .* 24 * 60 * 60 / (10^5))
    ex_n_arr = Interpolations.LinearInterpolation([t_i:t_e;], ex_n_arr .* 10^3)
    u_arr = Interpolations.LinearInterpolation([t_i:t_e;], fill(u, Int(t_e - t_i + 1)))

    return(u_arr, t_arr, irr_arr, ex_n_arr)
end

function solvekelp(t_i, nd, u, temp, irr, ex_n, lat, a_0, n_0, c_0, params="src/parameters/origional.jl", resp_model=1, dt=1)
    include(params)
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

    params = (u, temp, irr, ex_n, NormDeltaL, resp_model, dt)

    y_0 = vcat(a_0, n_0, c_0, 0)

    solver = ODEProblem(equations!, y_0, (t_i, t_i + nd), params)
    solution = solve(solver, RK4(), dt=dt, adaptive=false)#Please keep the RK4 algorithm otherwise the extreme caron limit needs to be changed

    results =
        DataFrame(area=[], nitrogen=[], carbon=[], gross_nitrate=[], time=[])
    for (ind, val) in enumerate(solution.u)
        push!(val, solution.t[ind])
        push!(results, (val))
    end

    return(solution, results)
end

end # module
