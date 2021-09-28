# Kelp.jl
[![DOI](https://zenodo.org/badge/383172934.svg)](https://zenodo.org/badge/latestdoi/383172934)[![Build Status](https://app.travis-ci.com/jagoosw/Kelp.jl.svg?branch=main)](https://app.travis-ci.com/jagoosw/Kelp.jl)[![codecov](https://codecov.io/gh/jagoosw/Kelp.jl/branch/main/graph/badge.svg?token=JG0D8UY2K8)](https://codecov.io/gh/jagoosw/Kelp.jl)

[Kelp.jl](https://github.com/jagoosw/Kelp.jl)  is an implimentation of the [Broch and Slagstad, 2012 model of the growth and composition of _Saccharina latissima_](https://link.springer.com/article/10.1007/s10811-011-9695-y).

The main way to solve a single frond is `Kelp.solvekelp` and grids can be solved by `Kelp.solvegrid`.

Changes from the stated parameter values in the paper are detailed in [changes.pdf](https://github.com/jagoosw/Kelp/blob/main/changes.pdf).

The package is not yet registered so to use, download this repository and then install the dependencies by executing (from this directory):
```
>julia
julia> import Pkg
julia> ] activate .
julia> instantiate
```
## Running a model
These examples for running the model both at a single point and on a grid can be found in [examples.ipynb](https://github.com/jagoosw/Kelp.jl/blob/main/examples/example.ipynb).

### Example 1 - Single Point


```julia
import Pkg; Pkg.activate("../")
using Kelp, Plots, Interpolations;pyplot();
```

Set initial conditions and parameters
- t_i - start day (days since january first) - needs to have the same date referance as the forcing data
- nd - number of days to run for
- lat - latitude as this effects the change in day length
- a_0,n_0,c_0 - initial area, nitrogen reserve and carbon reserve


```julia
t_i = 1.0;nd = 600;lat = 60
a_0 = 0.1;n_0 = 0.022;c_0 = 0.3;
```

Setting the forcing variables, these would normally be loaded from a data set but here will be generated
- time - day corresponding to forcing data in days since january first (year is arbitary). This **must** be a float rather than an integer or ODE solving library fails
- temp - temperature in degrees C
- no3 - nitrate concentration in mmol/m³
- irr - PAR irradiance in μmol photons/m^2/s
- u - water velocity in m/s


```julia
time = [1:2:800;]

temp = 6 * cos.((time .- 250) .* 2 .* pi ./ 365) .+ 8
no3 = (7 .* ((cos.(time .* 2 .* pi ./ 365) .+ 1) ./ 2).^3 .+ 0.1) ./ 1000
irr = 40 .* (sin.((time .+ 15) .* pi ./ 365).^10) .+ 1
u = repeat([0.15],400);
```


```julia
plot(layout=grid(1,3),size=(1000,250),legend=false)
plot!(time,temp,ylabel="Temp (°C)",sp=1)
plot!(time,no3,xlabel="Day",ylabel="Nitrate concentration (mmol/m³)",sp=2)
plot!(time,irr,ylabel="PAR (μmol photons/m²/s)",sp=3)
```




    
![png](output_6_0.png)
    



The forcing variables must be converted to interpolations for the kelp model to access them at arbitary time


```julia
temp_itp=Interpolations.LinearInterpolation(time,temp)
no3_itp=Interpolations.LinearInterpolation(time,no3)
irr_itp=Interpolations.LinearInterpolation(time,irr)
u_itp=Interpolations.LinearInterpolation(time,u);
```

Now the model can be run, the parameter file must be passed and in this run the resparation model proposed in Broch, 2013 is being used


```julia
solution, results = Kelp.solvekelp(t_i, nd, u_itp, temp_itp, irr_itp, no3_itp, lat, a_0, n_0, c_0, "../src/parameters/2013.jl",2);
```

Solutions contains the raw output of the ODE solver while results is refactored into a dataframe (this can optionally be turned off for an array to be returned)


```julia
results
```
*dataframe as a table that can't be rendered here*

It is useful to convert the results into total carbon and nitrogen masses (rather than the reserves that the model returns), this requires some of the parameters.


```julia
include("../src/parameters/2013.jl")
total_carbon = results.area .* K_A .* (results.carbon .+ C_struct)
total_nitrogen = results.area .* K_A .* (results.nitrogen .+ N_struct);
```


```julia
plot(layout=grid(1,3),size=(1000,250),legend=false)
plot!(results.time,results.area,ylabel="Area/dm²",sp=1)
plot!(results.time,total_carbon,xlabel="Day",ylabel="Total Carbon (gC)",sp=2)
plot!(results.time,total_nitrogen,ylabel="Total Nitrogen (gN)",sp=3)
```




    
![png](output_15_0.png)
    



### Example 2 - Grid

For a grid we must set initial conditions as with a single point


```julia
t_i = 1.0;nd = 300
a_0 = 0.1;n_0 = 0.022;c_0 = 0.3;
```

This time we will generate a 4d grid of input data for temp and no3 and 2d for nitrate. An aditional variable needs to be generated, either a 3d diffuse attenuation coefficient or 4d light attenuation. In this example the latter is used. Again this would usually be loaded from a file. It is benefitial to define all of these as constants as it drastically speeds up on larger grids.


```julia
const arr_lon=[45:50;]
const arr_lat=[55:65;]
const arr_dep=[0:5:75;]
const arr_t = [0:2:310;]

const arr_temp = permutedims(repeat(6 * cos.((arr_t .- 250) .* 2 .* pi ./ 365) .+ 8,1,length(arr_lon),length(arr_lat),length(arr_dep)),(2,3,4,1)).*permutedims(repeat(arr_lat./arr_lat[1],1,length(arr_lon),length(arr_dep),length(arr_t)),(2,1,3,4))
const arr_no3 = permutedims(repeat((7 .* ((cos.(arr_t .* 2 .* pi ./ 365) .+ 1) ./ 2).^3 .+ 0.1) ./ 1000,1,length(arr_lon),length(arr_lat),length(arr_dep)),(2,3,4,1)).*repeat(arr_lon./arr_lon[1],1,length(arr_lat),length(arr_dep),length(arr_t))
const arr_irr = permutedims(repeat(40 .* (sin.((arr_t .+ 15) .* pi ./ 365).^10) .+ 1,1,length(arr_lon),length(arr_lat)),(2,3,1)).*permutedims(repeat(arr_lat./arr_lat[1],1,length(arr_lon),length(arr_t)),(2,1,3))
const arr_beta = permutedims(repeat(reverse([0:1/(length(arr_dep)-1):1;]),1,length(arr_lon),length(arr_lat),length(arr_t)),(2,3,1,4))
const arr_u = permutedims(repeat([0.15],length(arr_t),length(arr_lon),length(arr_lat),length(arr_dep)),(2,3,4,1));
```

These grids are directly fed to the grid solver, which returns an array. irr can have its own time provided as satellite products often do not have the same temporal resolution as models. Additionally a fill value, in this case NaN, can be provided as they are often temporally sparse and need to be filtered.

This function automatically paralalises to however many threads you start julia with


```julia
@time results = Kelp.solvegrid(t_i, nd, a_0, n_0, c_0, arr_lon, arr_lat, arr_dep, arr_t, arr_no3, arr_temp, arr_u, (arr_irr, arr_t, NaN), (nothing, nothing, nothing), arr_beta, "../src/parameters/2013.jl", 2);
```
...

The output from this is an array with area/nitrogen reserve/carbon reserve/nitrate uptake in the first dimention, then lon,lat,dep,time in the others. We can extract the total carbon and nitrogen again:


```julia
total_carbon = results[1,:,:,:,:] .* K_A .* (results[3,:,:,:,:] .+ C_struct)
total_nitrogen = results[1,:,:,:,:] .* K_A .* (results[2,:,:,:,:] .+ N_struct);
```

We could plot these for a couple of points as a comparison for above:


```julia
plot(layout=grid(1,3),size=(1000,300))
results_time=[0:nd;]
for r=1:10
    i,j,k=rand([1:length(arr_lon);]),rand([1:length(arr_lat);]),rand([1:length(arr_dep);])
    plot!(results_time,results[1,i,j,k,:],sp=1,label="$(arr_lat[j])N, $(arr_lon[i])W, $(arr_dep[k])m")
    plot!(results_time,total_carbon[i,j,k,:],sp=2)
    plot!(results_time,total_nitrogen[i,j,k,:],sp=3)
end
plot!(ylabel="Area/dm²",sp=1,legend=:bottomright);plot!(sp=2,xlabel="Day",ylabel="Total Carbon (gC)",legend=false);plot!(sp=3,ylabel="Total Nitrogen (gN)",legend=false)
```




    
![png](output_26_0.png)
    



Or we could plot a heatmap of the surfaces:


```julia
hms=[
    heatmap(arr_lon,arr_lat,total_carbon[:,:,1,end]',color=cgrad(:bamako, rev=true),colorbar_title="Total Carbon (gC)"),
    heatmap(arr_lon,arr_lat,total_nitrogen[:,:,1,end]',color=:lajolla,colorbar_title="Total Nitrogen (gN)"),
    heatmap(arr_lon,arr_lat,total_carbon[:,:,1,end]'./total_nitrogen[:,:,1,end]',color=cgrad(:lapaz, rev=true),colorbar_title="Carbon:Nitrogen ratio")
]
plot!(hms...,layout=grid(1,3),size=(1000,200))
```

![png](output_28_0.png)

# Model Verification
The models outputs compared with figure 3 in Broch and Slagstad 2012 are shown below:
![B&S2012 Figure 3 equivalent.](img/paper_comparison.png)

The code to reproduce this can be found in the example folder.

The discrepancies may be down to the inaccuracy of reporting of the model parameters or the difficulty extracting the input data from the paper. 

Changes to the parameters from those published are detailed [here](https://github.com/jagoosw/Kelp/blob/master/changes.pdf).

## Check on the sidebar for function documentation