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

"""
    Kelp.equations!(y, params, t)

Equations outlined in the __Main Equations__ section of the paper to be solved by the ODE library.
- `y`: the current state of the system as a vector (area, nitrate reserve, carbon reserve).
- `params`: the variable parameters of the model:
    - `u_arr`: interpolation object of the water speed in time
    - `temp_arr`: interpolation object of the temperature in time
    - `irr_arr`: interpolation object of the irradianec in time
    - `ex_n_arr`: interpolation object of the external nitrate concentration in time
    - `NormDeltaL`: the normalised change in day length
    - `resp_model`: the choice of respiration model, 1 is the origional from the 2012 paper and 2 is the modified version in 2013
    - `dt`: the time step length, this is important as it is used in the "extreme carbon limit" part of the equations, see NB.
- `t`: the current time (with respect to the time in the interpolations)

Note:
These equations **must** be solved with an algorithm with fixed time steps and known, constant sub timestep lengths.
This is because the "extreme carbon limit" element is only implientable with known and fixed timesteps as the value of a 
must be changed to a particular value rather than changing the deriviative. This can only be done (within the framework of
the ODE library) by setting the derivitive to (X(next)-X(old))/dt
"""
function equations!(y::Vector{Float64}, params, t::Float64)
    a, n, c = y

    if a > 0
        if c <= C_min 
            a_check = a - a * (C_min - c) / C_struct
        else
            a_check = 0
        end

        u_arr, temp_arr, irr_arr, ex_n_arr, NormDeltaL, resp_model, dt = params

        u = u_arr(t)::Float64
        temp = temp_arr(t)::Float64
        irr = irr_arr(t)::Float64
        ex_n = ex_n_arr(t)::Float64

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
            @warn "Out of range temperature, $temp"
            f_temp = 0
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
        da, dn, dc, j = 0, 0, 0, 0 
        @warn "Area reached 0"
    end
    return (vcat(da, dn, dc, j * a))
end

"""
    Kelp.defaults(t_i, t_e, u)
Generates "default" or anayltical values for the water speed, temperature, irradiance and external nitrogen for testing.

Parameters:
- `t_i`: start time
- `t_e`: end time
- `u`: your chosen water speed
"""
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

"""
    Kelp.solvekelp(t_i, nd, u, temp, irr, ex_n, lat, a_0, n_0, c_0, params="src/parameters/origional.jl", resp_model=1, dt=1)
Solves the model for some set of parametetrs and returns the ODE library solution as well as a dataframe of the useful results.

Parameters:
- `t_i`: th estart time (in days since the start of the interpolation objects "day zero")
- `nd`: the number of days to run for
- `u`: interpolation object (IO) of water speed
- `temp`: IO of temperature
- `ex_n`: IO of external nitrate concentration
- `lat`: latitude, relivant for the change of day length
- `a_0`: initial area
- `n_0`: initial nitrogen reserve (gN/gSW)
- `c_0`: initial carbon reserve (gC/gSW)
- `params`: string of the path to a parameters file, defaults to the 2012 values. Also supplied is 2013 in src/parameters/2013.jl or you can copy and vary them
- `resp_model`: choice of respiration model, 1 (default) uses the 2012 version and 2 uses the modifcations from the 2013 paper
- `dt`: the time step size to use (see equations! note), default is 1 day (seems small enough)

Returns:
- `solution`: the ODE library solution
- `results`: dataframe of area/nitrogen reserve/carbon reserve/total nitrate update. All others useful quantities can be easily derived.
"""
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
    solution = solve(solver, RK4(), dt=dt, adaptive=false)# Please keep the RK4 algorithm otherwise the extreme caron limit needs to be changed

    results =
        DataFrame(area=[], nitrogen=[], carbon=[], gross_nitrate=[], time=[])
    for (ind, val) in enumerate(solution.u)
        push!(val, solution.t[ind])
        push!(results, (val))
    end

    return(solution, results)
end

"""
    Kelp.solvegrid(t_i, nd, a_0, n_0, c_0, arr_lon, arr_lat, arr_dep, arr_time, no3, temp, u, par_data, kd_data, params="src/parameters/origional.jl", resp_model=1, dt=1)
Solve the model for a (spacially) fixed grid of inputs.

Parameters:
- `t_i`: start day (w.r.t. t=0 in time)
- `nd`: number of days to run for
- `a_0`: initial area/dm^2
- `n_0`: initial nitrogen reserve/gN/gSW
- `c_0`: initial carbon reserve/gC/gSW
- `arr_lon`: longitudes for no3 and temp
- `arr_lat`: latitudes for no3 and temp
- `arr_dep`: depths for no3 and temp
- `arr_time`: time for no3 and temp
- `no3`: array of no3 concentration in lon,lat,depth,time
- `temp`: array of temp concentration in lon,lat,depth,time
- `u`: array of water speed in lon,lat,depth,time
- `par_data`: array of:
    - `par values in lon,lat,time
    - `corresponding time
    - `par fill value
- `kd_data`: array of:
    - kd values in lon,lat,time
    - corresponding time
    - kd fill value
- `params`: string of the path to a parameters file, defaults to the 2012 values. Also supplied is 2013 in src/parameters/2013.jl or you can copy and vary them
- `resp_model`: choice of respiration model, 1 (default) uses the 2012 version and 2 uses the modifcations from the 2013 paper
- `dt`: the time step size to use (see equations! note), default is 1 day (seems small enough)

Returns: results as an array of (area/nitrogen/carbon/nitrate update, lon, lat, depth, time)

Notes:
par and kd need their own time coordinates and fill value because they come from satelite observation which are 
temporally sparse so need to be checked and interpolated in time for each point. On the other hand temp and no3
(that I'm using) are from Copurnicus' models so if a point has a value at some time it will at all times.

no3,temp and u need to be of the same shape and size and with the values corresponding to the same position/time.
"""
function solvegrid(t_i::Float64, nd::Int, a_0::Float64, n_0::Float64, c_0::Float64, arr_lon, arr_lat, arr_dep, arr_time, no3::Array{Float64,4}, temp::Array{Float64,4}, u::Array{Float64,4}, par_data, kd_data, params::String="src/parameters/origional.jl", resp_model::Int=1, dt=1)
    # Would like to annotate type for the others but for some reason they making tuples of "Number" doesn't isn't satisfied
    par, par_t, par_fill = par_data;kd, kd_t, kd_fill = kd_data

    all_results::Array{Float64,5} = repeat([NaN], 4, length(arr_lon), length(arr_lat), length(arr_dep), floor(Int, nd) + 1);
    points::Array{Vector{Int64},3} = collect.(Iterators.product([1:length(arr_lon);], [1:length(arr_lat);], [1:length(arr_dep);]));

    Threads.@threads for (i, j, k) in points
        lat = arr_lat[j];depth = arr_dep[k]

        no3_vals = no3[i,j,k,:]
        temp_vals = temp[i,j,k,:]
        u_vals = u[i,j,k,:]
        if (!isnan(no3_vals[1])) & (!isnan(temp_vals[1])) & (!isnan(u_vals[1])) & all(n > 0 for n in no3_vals)
            no3_itp = Interpolations.LinearInterpolation(arr_time, no3_vals)
            temp_itp = Interpolations.LinearInterpolation(arr_time, temp_vals)
            u_itp = Interpolations.LinearInterpolation(arr_time, u_vals)
            
            par_vals_raw = par[i,j,:]
            kd_vals_raw = kd[i,j,:]

            par_vals, par_t_vals = extract_valid(par_vals_raw, par_t, par_fill)
            kd_vals, kd_t_vals = extract_valid(kd_vals_raw, kd_t, kd_fill)

            if (length(par_vals) > 6) & (length(kd_vals) > 6)
                par_itp = Interpolations.LinearInterpolation(par_t_vals, par_vals, extrapolation_bc=Flat())
                kd_itp = Interpolations.LinearInterpolation(kd_t_vals, kd_vals, extrapolation_bc=Flat())
                irr_itp = Interpolations.LinearInterpolation(arr_time, par_itp.(arr_time) .* exp.(-kd_itp.(arr_time) .* depth), extrapolation_bc=Flat())

                solution, results = Kelp.solvekelp(t_i, nd, u_itp, temp_itp, irr_itp, no3_itp, lat, a_0, n_0, c_0, params, resp_model, dt);
                
                all_results[1,i,j,k,:] = results.area;
                all_results[2,i,j,k,:] = results.nitrogen;
                all_results[3,i,j,k,:] = results.carbon;
                all_results[4,i,j,k,:] = results.gross_nitrate
            end
        end
    end
    return all_results
end

"""
    Kelp.extract_valid(raw,raw_time,fill)
Extracts the valid values from an array by checking against a fill value and returns the valids and corresponding time.

Parameters:
- `raw`: the array to check
- `raw_time`: the corresponding time array
- `fill`: the fill value to check against

Returns:
- `vals`: the filtered values
- `time`: corresponding times
"""
function extract_valid(raw, raw_time, fill)
    vals, times = [], []
    for (ind, val) in enumerate(raw)
        if val != fill
            push!(vals, val)
            push!(times, raw_time[ind])
        end
    end
    return (vals, times)
end
end # module
"""
    get_int(val, list, tol)
Function that finds the index in the list with the closest value to val. Error is thrown if no result is within tollerance, tol.

Parameters:
- `val`: the value searching for
- `list`: the list to search
- `tol`: tollerance of search

Returns:
- `index of closest valuea
"""
function get_ind(val, list, tol)
    result = findmin(abs.(list .- val))
    ind = result[2]
    if result[1] > tol
        closest = list[result[2]]
        throw("No indicie could be found within tollerance for the requested value, requested was $val, closes was $closest at $ind")
    else
        return ind
    end
end

"""
    interp_deps(arr, origional_depths, desired_depths, invalid_val)
Function to interpolate a 4D array in the last dimension. Useful for lineaising the depth steps of the arrays
from Copurnicus as they have increasingly corse step sizes.

Parameters:
- `arr`: the array to interpolate where the 3rd dimension is the one to interpolate along
- `origional_depths`: the values of the 3rd dimensions coordinates (the origional depths of the data)
- `desired_depths`: the desired coordinates/depths
- `invalid_val`: the fill value to replace with NaN in the array (as we are searching it anyway)

Returns: new array with interpolated 3rd dimension
"""
function interp_deps(arr, origional_depths, desired_depths, invalid_val)
    arr_size = size(arr)
    new_arr = repeat([NaN], arr_size[1], arr_size[2], length(desired_depths), arr_size[4]);

    points = collect.(Iterators.product([1:arr_size[1];], [1:arr_size[2];], [1:arr_size[4];]));
    Threads.@threads for (i, j, k) in points
        dep_arr = arr[i,j,:,k]
        usable = [];
        finished = false
        for (ind, val) in enumerate(dep_arr)
            if (finished == false) & (!isnan(val)) & (round(val) != round(invalid_val))
                push!(usable, val)
            else
                finished = true
            end
        end
        if length(usable) > 1
            deps = origional_depths[1:length(usable)]
            dep_itp = Interpolations.LinearInterpolation(deps, usable, extrapolation_bc=Flat())

            extent = findmin(abs.(desired_depths .- deps[end]))[2]
            if extent != length(desired_depths)
                search_deps = desired_depths[1:extent - 1]
            else
                search_deps = desired_depths
            end
            new_arr[i,j,1:length(search_deps),k] = dep_itp.(search_deps)
        end
    end

    return new_arr
end