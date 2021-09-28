# Example 1 - Single Point


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




<div class="data-frame"><p>601 rows × 5 columns</p><table class="data-frame"><thead><tr><th></th><th>area</th><th>nitrogen</th><th>carbon</th><th>gross_nitrate</th><th>time</th></tr><tr><th></th><th title="Any">Any</th><th title="Any">Any</th><th title="Any">Any</th><th title="Any">Any</th><th title="Any">Any</th></tr></thead><tbody><tr><th>1</th><td>0.1</td><td>0.022</td><td>0.3</td><td>0.0</td><td>1.0</td></tr><tr><th>2</th><td>0.109134</td><td>0.0189379</td><td>0.267903</td><td>8.52351e-8</td><td>2.0</td></tr><tr><th>3</th><td>0.116672</td><td>0.0167756</td><td>0.247928</td><td>3.68739e-7</td><td>3.0</td></tr><tr><th>4</th><td>0.122468</td><td>0.0152978</td><td>0.237388</td><td>8.09489e-7</td><td>4.0</td></tr><tr><th>5</th><td>0.126655</td><td>0.014318</td><td>0.23388</td><td>1.36651e-6</td><td>5.0</td></tr><tr><th>6</th><td>0.129536</td><td>0.0136848</td><td>0.235377</td><td>2.00431e-6</td><td>6.0</td></tr><tr><th>7</th><td>0.131451</td><td>0.0132832</td><td>0.240287</td><td>2.69541e-6</td><td>7.0</td></tr><tr><th>8</th><td>0.132698</td><td>0.0130318</td><td>0.247446</td><td>3.42009e-6</td><td>8.0</td></tr><tr><th>9</th><td>0.133506</td><td>0.0128758</td><td>0.256043</td><td>4.16497e-6</td><td>9.0</td></tr><tr><th>10</th><td>0.134031</td><td>0.0127794</td><td>0.265536</td><td>4.92123e-6</td><td>10.0</td></tr><tr><th>11</th><td>0.134378</td><td>0.01272</td><td>0.27557</td><td>5.68316e-6</td><td>11.0</td></tr><tr><th>12</th><td>0.134615</td><td>0.0126834</td><td>0.28592</td><td>6.4471e-6</td><td>12.0</td></tr><tr><th>13</th><td>0.134782</td><td>0.0126608</td><td>0.296442</td><td>7.21063e-6</td><td>13.0</td></tr><tr><th>14</th><td>0.134907</td><td>0.0126468</td><td>0.307046</td><td>7.97216e-6</td><td>14.0</td></tr><tr><th>15</th><td>0.135005</td><td>0.0126382</td><td>0.317676</td><td>8.73057e-6</td><td>15.0</td></tr><tr><th>16</th><td>0.135087</td><td>0.0126328</td><td>0.328295</td><td>9.48508e-6</td><td>16.0</td></tr><tr><th>17</th><td>0.135158</td><td>0.0126295</td><td>0.338884</td><td>1.02351e-5</td><td>17.0</td></tr><tr><th>18</th><td>0.135223</td><td>0.0126274</td><td>0.349427</td><td>1.09801e-5</td><td>18.0</td></tr><tr><th>19</th><td>0.135284</td><td>0.0126261</td><td>0.359918</td><td>1.17197e-5</td><td>19.0</td></tr><tr><th>20</th><td>0.135342</td><td>0.0126253</td><td>0.370352</td><td>1.24535e-5</td><td>20.0</td></tr><tr><th>21</th><td>0.135398</td><td>0.0126248</td><td>0.380726</td><td>1.31812e-5</td><td>21.0</td></tr><tr><th>22</th><td>0.135452</td><td>0.0126244</td><td>0.391039</td><td>1.39026e-5</td><td>22.0</td></tr><tr><th>23</th><td>0.135506</td><td>0.0126242</td><td>0.401291</td><td>1.46173e-5</td><td>23.0</td></tr><tr><th>24</th><td>0.135559</td><td>0.012624</td><td>0.411481</td><td>1.5325e-5</td><td>24.0</td></tr><tr><th>25</th><td>0.13561</td><td>0.0126239</td><td>0.42161</td><td>1.60256e-5</td><td>25.0</td></tr><tr><th>26</th><td>0.135662</td><td>0.0126238</td><td>0.431679</td><td>1.67187e-5</td><td>26.0</td></tr><tr><th>27</th><td>0.135712</td><td>0.0126237</td><td>0.441689</td><td>1.74041e-5</td><td>27.0</td></tr><tr><th>28</th><td>0.135762</td><td>0.0126237</td><td>0.451641</td><td>1.80815e-5</td><td>28.0</td></tr><tr><th>29</th><td>0.135812</td><td>0.0126236</td><td>0.461535</td><td>1.87509e-5</td><td>29.0</td></tr><tr><th>30</th><td>0.135861</td><td>0.0126235</td><td>0.471374</td><td>1.94119e-5</td><td>30.0</td></tr><tr><th>&vellip;</th><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td></tr></tbody></table></div>



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
    



# Example 2 - Grid

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

    ┌ Info: At level 1
    └ @ Kelp /home/jago/Documents/Projects/Kelp.jl/src/Kelp.jl:337
    ┌ Info: At level 2
    └ @ Kelp /home/jago/Documents/Projects/Kelp.jl/src/Kelp.jl:337
    ┌ Info: At level 3
    └ @ Kelp /home/jago/Documents/Projects/Kelp.jl/src/Kelp.jl:337
    ┌ Info: At level 4
    └ @ Kelp /home/jago/Documents/Projects/Kelp.jl/src/Kelp.jl:337
    ┌ Info: At level 5
    └ @ Kelp /home/jago/Documents/Projects/Kelp.jl/src/Kelp.jl:337
    ┌ Info: At level 6
    └ @ Kelp /home/jago/Documents/Projects/Kelp.jl/src/Kelp.jl:337
    ┌ Info: At level 7
    └ @ Kelp /home/jago/Documents/Projects/Kelp.jl/src/Kelp.jl:337
    ┌ Info: At level 8
    └ @ Kelp /home/jago/Documents/Projects/Kelp.jl/src/Kelp.jl:337
    ┌ Info: At level 9
    └ @ Kelp /home/jago/Documents/Projects/Kelp.jl/src/Kelp.jl:337
    ┌ Info: At level 10
    └ @ Kelp /home/jago/Documents/Projects/Kelp.jl/src/Kelp.jl:337
    ┌ Info: At level 11
    └ @ Kelp /home/jago/Documents/Projects/Kelp.jl/src/Kelp.jl:337
    ┌ Info: At level 12
    └ @ Kelp /home/jago/Documents/Projects/Kelp.jl/src/Kelp.jl:337
    ┌ Info: At level 13
    └ @ Kelp /home/jago/Documents/Projects/Kelp.jl/src/Kelp.jl:337
    ┌ Info: At level 14
    └ @ Kelp /home/jago/Documents/Projects/Kelp.jl/src/Kelp.jl:337
    ┌ Info: At level 15
    └ @ Kelp /home/jago/Documents/Projects/Kelp.jl/src/Kelp.jl:337
    ┌ Info: At level 16
    └ @ Kelp /home/jago/Documents/Projects/Kelp.jl/src/Kelp.jl:337


     32.622805 seconds (135.82 M allocations: 4.460 GiB, 2.42% gc time, 25.27% compilation time)


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
    




```julia

```
