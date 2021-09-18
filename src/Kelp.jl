module Kelp
using RecursiveArrayTools, DiffEqBase, OrdinaryDiffEq, DiffEqCallbacks, Roots, Interpolations, DataFrames
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
    Kelp.eval_μ(a,n,c,temp,λ)

Solves Equation 2, the specific growth rate.

Parameters:
- `a`: area /dm^2
- `n`: nitrate reserve /gN/gSW
- `c`: carbon reserve /gC/gSW
- `temp`: temperature /°C
- `λ`: normalised change in day length

Returns: specific growth rate

Notes:
- Names eval_μ because it looks better in Kelp.equations to have μ and you can't reuse it
"""
eval_μ(a,n,c,temp,λ) = f_area(a) * f_temp(temp) * f_photo(λ) * min(1 - N_min / n, 1 - C_min / c)

"""
    Kelp.f_area(a)

Solves Equation 3, the effect of area on growth

Parameters: `a`, area /dm^2

Returns: effect of area on growth
"""
f_area(a) = m_1 * exp(-(a / A_0)^2) + m_2

"""
    Kelp.f_temp(temp)

Solves Equation 4, the effect of temperature on growth

Parameters: `temp`: temperature /°C

Returns: effect of temperature on growth
"""
function f_temp(temp)
    if -1.8 <= temp < 10 
        return 0.08 * temp + 0.2
    elseif 10 <= temp <= 15
        return 1
    elseif 15 < temp <= 19
        return 19 / 4 - temp / 4
    elseif temp > 19
        return 0
    else
        return 0
    end
end

"""
    Kelp.f_photo(λ)

Solves Equation 5, the seasonal influence on growth rate

Parameters: `λ`: normalised change in day length

Returns: seasonal influence on growth
"""
f_photo(λ) = a_1 * (1 + sign(λ) * abs(λ)^.5) + a_2

"""
    Kelp.ν(a)

Solves Equation 6, the specific frond erosion rate.

Parameters: `a`: area /dm^2

Returns: specific frond erosion rate
"""
ν(a) = 1e-6 * exp(ϵ * a) / (1 + 1e-6 * (exp(ϵ * a) - 1))

"""
    Kelp.eval_j(ex_n,n,u)

Solves Equation 8, the specific nitrate uptake rate.

Parameters:
- `ex_n`: external nitrate concentration /mmol/m^3
- `n`: nitrogen reserve /gN/gSW
- `u`: water speed /m/s

Returns: specific nitrate uptake rate

Notes:
- Names eval_j because it looks better in Kelp.equations to have j and you can't reuse it
"""
eval_j(ex_n,n,u) = J_max * (ex_n / (K_X + ex_n)) * ((N_max - n) / (N_max - N_min)) * (1 - exp(-u / U_0p65))

"""
    Kelp.p(temp,irr)

Solves Equation 10, the gross photosynthesis.

Parameters:
- `temp`: temperature /°C
- `irr`: irradiance /μmol photons/m^2/s

Returns: gross photosynthesis function
"""
function p(temp, irr)
    p_max = P_1 * exp(T_AP / T_P1 - T_AP / (temp + 273.15)) / (1 + exp(T_APL / (temp + 273.15) - T_APL / T_PL) + exp(T_APH / T_PH - T_APH / (temp + 273.15)))
    β_func(x) = p_max - (α * I_sat / log(1 + α / x)) * (α / (α + x)) * (x / (α + x))^(x / α)
    β = find_zero(β_func, (0, 0.1), Bisection())
    p_s = α * I_sat / log(1 + α / β)
    return p_s * (1 - exp(-α * irr / p_s)) * exp(-β * irr / p_s) 
end

"""
    Kelp.r(temp,μ,j,resp_model)

Solves Equation 14 (2012) if resp_model=1, or Equation 2 (2013) if resp_model=2

Parameters:
- `temp`: temperature /°C
- `μ`: specific growth rate
- `j`: specific nitrate uptake rate
- `resp_model`: choice of resparation model (see description)

Returns: respiration function
"""
function r(temp, μ, j, resp_model)
    if resp_model == 1
        return R_1 * exp(T_AR / T_R1 - T_AR / (temp + 273.15)) # temperature dependent respiration
    elseif resp_model == 2
        return (R_A * (μ / μ_max + j / J_max) + R_B) * exp(T_AR / T_R1 - T_AR / (temp + 273.15))
    end
end

e(c) = 1 - exp(γ * (C_min - c))

"""
    Kelp.equations(t, a, n, c, u, temp, irr, ex_n, λ, resp_model, dt)

Solves the papers main equations.

Parameters:
- `a`: area /dm^2
- `n`: nitrate reserve /gN/gSW
- `c`: carbon reserve /gC/gSW
- `u`: water speed /m/s
- `temp`: temperature /°C
- `irr`: irradiance / μmol photons/m^2/s 
- `ex_n`: external nitrate concentration /mmol m^3
- `λ`: normalised change in day length
- `resp_model`: choice of resparation model (see Kelp.r)
- `dt`: timestep length /days

Returns:
- `da`: results of equation 1, the area rate
- `dn`: results of equation 7, the nitrogen rate
- `dc`: results of equation 9, the carbon rate
- `j`: the specific nitrate uptake rate
"""
function equations(a, n, c, u, temp, irr, ex_n, λ, resp_model, dt)
    μ = eval_μ(a, n, c, temp, λ)
    j = eval_j(ex_n, n, u)

    # Equation 1
    da = (μ - ν(a)) * a
    # Equation 7
    dn = j / K_A - μ * (n + N_struct)
    # Equation 9
    dc = (p(temp, irr) * (1 - e(c)) - r(temp, μ, j, resp_model)) / K_A - μ * (c + C_struct)

    return da, dn, dc, j
end

function extremecarbon(y, t, integrator)
    if y[3]<C_min
        y[1]-=y[1]*(C_min-y[3])/C_struct
        y[3]=C_min
        println("ext carb")
    end
end
"""
    Kelp.solver!(y, params, t)

Interface for OrdinaryDiffEq library, extracts current enviromental variable values and solves equations for each timestep.

Parameters:
- `y`: the current state of the system as a vector (area, nitrate reserve, carbon reserve).
- `params`: the variable parameters of the model:
    - `u_arr`: interpolation object of the water speed in time
    - `temp_arr`: interpolation object of the temperature in time
    - `irr_arr`: interpolation object of the irradianec in time
    - `ex_n_arr`: interpolation object of the external nitrate concentration in time
    - `λ_arr`: array of the normalised change in day length
    - `resp_model`: the choice of respiration model, 1 is the origional from the 2012 paper and 2 is the modified version in 2013
    - `dt`: the time step length, this is important as it is used in the "extreme carbon limit" part of the equations, see NB.
- `t`: the current time (with respect to the time in the interpolations)

Returns: array of da,dn,dc,j*a - the rate of a/n/c and the nitrate uptake rate

Note:
These equations **must** be solved with an algorithm with fixed time steps and known, constant sub timestep lengths.
This is because the "extreme carbon limit" element is only implientable with known and fixed timesteps as the value of a 
must be changed to a particular value rather than changing the deriviative. This can only be done (within the framework of
the ODE library) by setting the derivitive to (X(next)-X(old))/dt
"""
function solver!(y::Vector{Float64}, params, t::Float64)
    a, n, c = y

    if a > 0
        u_arr, temp_arr, irr_arr, ex_n_arr, λ_arr, resp_model, dt = params

        u = u_arr(t)::Float64
        temp = temp_arr(t)::Float64
        irr = irr_arr(t)::Float64
        ex_n = ex_n_arr(t)::Float64

        d = trunc(Int, mod(floor(t), 365) + 1) 
        λ = λ_arr[d] 

        da, dn, dc, j = equations(a, n, c, u, temp, irr, ex_n, λ, resp_model, dt)

    else
        da, dn, dc, j = 0, 0, 0, 0 
    end
    return (vcat(da, dn, dc, j * a))
end

"""
    Kelp.solvekelp(t_i, nd, u, temp, irr, ex_n, lat, a_0, n_0, c_0, params="src/parameters/origional.jl", resp_model=1, dt=1, dataframe=true)

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
- `dataframe`: output as a dataframe, default to true. Alternative is an array (faster)
- `λ_arr`: array of normalised change of day length, defaults to nothing which generates the default one

Returns:
- `solution`: the ODE library solution
- `results`: dataframe or array of area/nitrogen reserve/carbon reserve/total nitrate update. All others useful quantities can be easily derived.
"""
function solvekelp(t_i, nd, u, temp, irr, ex_n, lat, a_0, n_0, c_0, params="../src/parameters/origional.jl", resp_model=1, dt=1, dataframe=true, λ_arr=nothing)
    include(params)
    if λ_arr==nothing
        λ_arr = gen_λ(lat)
    end
    params = (u, temp, irr, ex_n, λ_arr, resp_model, dt)

    y_0 = vcat(a_0, n_0, c_0, 0)

    solver = ODEProblem(solver!, y_0, (t_i, t_i + nd), params)
    solution = OrdinaryDiffEq.solve(solver, Vern8(),callback=FunctionCallingCallback(extremecarbon,func_everystep=true))# Please keep the RK4 algorithm otherwise the extreme caron limit needs to be changed

    if dataframe == true
        results = DataFrame(area=[], nitrogen=[], carbon=[], gross_nitrate=[], time=[])
        for (ind, val) in enumerate(solution.u)
            push!(val, solution.t[ind])
            push!(results, (val))
        end
        return(solution, results)
    else
        println(solution.t)
        return vcat(transpose.(solution.u)...)
    end
end

"""
    Kelp.solvegrid(t_i, nd, a_0, n_0, c_0, arr_lon, arr_lat, arr_dep, arr_time, no3, temp, u, par_data, kd_data, att = nothing, params="src/parameters/origional.jl", resp_model=1, dt=1, progress=true)
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
- `att`: array of light attenuation coefficients (PAR(z)=PAR(z=0)*att) in lon,lat,depth,time. Defaults to nothing and kd is used instead
- `params`: string of the path to a parameters file, defaults to the 2012 values. Also supplied is 2013 in src/parameters/2013.jl or you can copy and vary them
- `resp_model`: choice of respiration model, 1 (default) uses the 2012 version and 2 uses the modifcations from the 2013 paper
- `dt`: the time step size to use (see equations! note), default is 1 day (seems small enough)
- `progress`: option to update progress at each level (when multithreading not accurate but useful), defaults to true

Returns: results as an array of (area/nitrogen/carbon/nitrate update, lon, lat, depth, time)

Notes:
par and kd need their own time coordinates and fill value because they come from satelite observation which are 
temporally sparse so need to be checked and interpolated in time for each point. On the other hand temp and no3
(that I'm using) are from Copurnicus' models so if a point has a value at some time it will at all times.

no3,temp and u need to be of the same shape and size and with the values corresponding to the same position/time.
"""
function solvegrid(t_i::Float64, nd::Int, a_0::Float64, n_0::Float64, c_0::Float64, arr_lon, arr_lat, arr_dep, arr_time, no3::Array{Float64,4}, temp::Array{Float64,4}, u::Array{Float64,4}, par_data, kd_data, att=nothing, params::String="../src/parameters/origional.jl", resp_model::Int=1, dt=1, progress::Bool=true)
    # Would like to annotate type for the others but for some reason they making tuples of "Number" doesn't isn't satisfied
    par, par_t, par_fill = par_data;kd, kd_t, kd_fill = kd_data

    all_results::Array{Float64,5} = repeat([NaN], 4, length(arr_lon), length(arr_lat), length(arr_dep), round(Int, nd / dt + 1));
    points::Array{Vector{Int64},3} = collect.(Iterators.product([1:length(arr_lon);], [1:length(arr_lat);], [1:length(arr_dep);]));

    Threads.@threads for (i, j, k) in points
        if (progress == true) & (i == 1) & (j == 1)
            @info("At level $k")
        end

        lat = arr_lat[j];depth = arr_dep[k]

        no3_vals = no3[i,j,k,:]
        temp_vals = temp[i,j,k,:]
        u_vals = u[i,j,k,:]
        if (!isnan(no3_vals[1])) & (!isnan(temp_vals[1])) & (!isnan(u_vals[1])) & all(n > 0 for n in no3_vals)
            no3_itp = Interpolations.LinearInterpolation(arr_time, no3_vals)
            temp_itp = Interpolations.LinearInterpolation(arr_time, temp_vals)
            u_itp = Interpolations.LinearInterpolation(arr_time, u_vals)
            
            par_vals_raw = par[i,j,:]
            par_vals, par_t_vals = extract_valid(par_vals_raw, par_t, par_fill)
            
            if att == nothing
                kd_vals_raw = kd[i,j,:]
                kd_vals, kd_t_vals = extract_valid(kd_vals_raw, kd_t, kd_fill)
            else
                att_vals = att[i,j,k,:]
            end

            if (length(par_vals) > 6)# & (att!=nothing) | (length(kd_vals) > 6)
                par_itp = Interpolations.LinearInterpolation(par_t_vals, par_vals, extrapolation_bc=Flat())
                if att == nothing
                    kd_itp = Interpolations.LinearInterpolation(kd_t_vals, kd_vals, extrapolation_bc=Flat())
                    irr_itp = Interpolations.LinearInterpolation(arr_time, par_itp.(arr_time) .* exp.(-kd_itp.(arr_time) .* depth), extrapolation_bc=Flat())
                else
                    irr_itp = Interpolations.LinearInterpolation(arr_time, par_itp.(arr_time) .* att_vals, extrapolation_bc=Flat())
                end
                solution = Kelp.solvekelp(t_i, nd, u_itp, temp_itp, irr_itp, no3_itp, lat, a_0, n_0, c_0, params, resp_model, dt, false);

                all_results[:,i,j,k,1:size(solution)[1]] = reshape(solution', size(solution)[2], 1, 1, 1, size(solution)[1])
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
        if (val != fill) & (!isnan(val))
            push!(vals, val)
            push!(times, raw_time[ind])
        end
    end
    return (vals, times)
end

"""
    Kelp.gen_λ(lat)

Generates λ as described in Model desctiptions/Main equations/Photoperiodic effect (page 763) in the origional paper.

Parameter: `lat`, the latitude

Returns: array of normalised change in day length
"""
function gen_λ(lat)
    θ = 0.2163108 .+ 2 .* atan.(0.9671396 .* tan.(0.00860 .* ([1:365;] .- 186)))
    δ = asin.(0.39795 .* cos.(θ))
    p = 0.8333
    d = 24 .- (24 / pi) .* acos.((sin(p * pi / 180) .+ sin(lat * pi / 180) .* sin.(δ)) ./ (cos(lat * pi / 180) .* cos.(δ)))
    λ = diff(d)
    push!(λ, λ[end])
    return λ ./ findmax(λ)[1]
end
"""
    Kelp.get_ind(val, list, tol)
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
    Kelp.interp_deps(arr, origional_depths, desired_depths, invalid_val)
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
end # module
