using Plots, Dates, Measures, Interpolations, RollingFunctions

include("../src/Kelp.jl")

#collected from fig2 of origional paper with https://www.graphreader.com/
temp_t=[20.956,31.796,43.358,65.036,67.927,82.38,86.715,99.723,111.285,124.292,147.416,209.562,238.467,296.277,328.073,356.978,394.555,412.62]
temp=[12.255,14.387,13.177,12.947,12.082,12.313,11.333,11.103,9.663,10.181,6.609,3.383,2.461,8.107,12.658,15.251,14.617,13.811]
no3_t=[26.737,79.489,128.628,203.781,228.35,258.701,268.818,299.168,325.182,399.613]
no3=[0.198,0.428,5.136,7.045,2.798,5.333,0.23,0.428,0.132,0.132]
irr_t=[22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,342,343,344,345,346,347,348,349,350,351,352,353,354,355,356,357,358,359,360,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375,376,377,378,379,380,381,382,383,384,385,386,387,388,389,390,391,392,393,394,395,396]
irr=[8.034,7.878,7.639,7.399,7.16,6.92,8.948,11.36,13.773,13.64,12.202,10.765,9.328,7.89,6.453,5.466,5.671,5.877,6.082,6.287,6.493,6.698,6.903,6.736,5.658,4.58,3.622,3.787,3.951,4.115,4.28,4.444,4.608,4.772,4.86,4.347,3.833,3.32,2.807,2.486,2.568,2.65,2.732,2.814,2.896,2.978,3.06,3.028,2.771,2.514,2.258,2.001,1.744,1.488,1.475,1.475,1.475,1.475,1.475,1.475,1.475,1.475,1.475,1.571,1.827,2.084,2.259,2.053,1.848,1.643,1.437,1.232,1.027,0.821,0.656,0.656,0.656,0.656,0.656,0.656,0.656,0.656,0.656,0.656,0.656,0.656,0.656,0.656,0.656,0.656,0.654,0.644,0.635,0.626,0.616,0.607,0.598,0.588,0.579,0.57,0.56,0.551,0.542,0.532,0.523,0.514,0.504,0.495,0.483,0.47,0.458,0.445,0.432,0.419,0.406,0.393,0.381,0.368,0.355,0.342,0.329,0.343,0.36,0.377,0.394,0.411,0.429,0.446,0.463,0.48,0.497,0.514,0.531,0.548,0.565,0.583,0.6,0.617,0.634,0.651,0.693,0.744,0.796,0.847,0.898,0.95,1.001,1.052,1.104,1.155,1.206,1.258,1.309,1.465,1.626,1.788,1.949,2.11,2.272,2.433,2.594,2.756,2.917,3.078,3.258,3.443,3.627,3.812,3.997,4.182,4.367,4.551,3.982,3.212,2.442,2.254,2.459,2.664,2.87,3.075,3.28,3.486,3.691,3.896,4.102,4.349,4.595,4.842,5.088,5.334,5.581,5.827,6.069,6.172,6.274,6.377,6.48,6.582,6.685,7,7.431,7.862,8.293,8.724,9.155,9.587,10.018,9.824,9.311,8.798,8.284,7.771,7.258,6.824,6.602,6.379,6.157,5.934,5.712,5.49,5.267,5.045,4.822,5.003,5.363,5.722,6.081,6.441,6.8,7.039,7.004,6.97,6.936,6.902,6.868,6.833,6.799,6.765,6.731,6.696,6.662,6.628,6.594,6.56,7.662,8.843,10.023,10.633,11.043,11.454,11.865,12.275,12.686,12.805,12.394,11.984,11.573,11.162,10.752,10.341,12.53,14.806,17.082,19.358,21.633,23.909,26.185,28.461,30.736,32.303,32.748,33.193,33.638,34.083,34.528,34.973,35.418,35.862,36.307,35.759,34.972,34.184,33.397,32.636,33.423,34.21,34.997,35.784,36.567,37.337,38.107,38.752,35.631,32.51,29.389,26.268,23.147,20.026,16.904,14.156,18.742,23.328,27.914,32.499,34.524,29.185,23.847,19.117,23.224,27.331,31.437,35.544,36.343,28.386,20.43,13.143,13.492,13.841,14.19,14.539,14.888,15.237,15.586,16.001,17.028,18.054,19.081,21.035,23.088,25.142,27.195,28.742,26.312,23.883,21.453,19.023,17.024,15.929,14.834,13.738,12.643,12.281,12.564,12.846,13.128,13.411,13.693,13.855,13.307,12.76,12.212,11.665,11.226,10.987,10.747,10.508,10.268,10.029,9.789,9.55,9.31,9.07,8.838,8.607,8.376,8.145,7.914,7.683,7.511,7.434,7.357,7.28,7.203,7.126,7.049]
no3.*=0.001#I don't think the units are correctly recorded in the paper
irr.*=10^6/(24*60*60)
temp_t.+=7*30;no3_t.+=7*30;irr_t.+=7*30

t_i=8.0*30
nd=365
lat=60.25
u=0.15

u_arr = Interpolations.LinearInterpolation([t_i:t_i+nd;], fill(u, Int(nd + 1)))

temp_arr=Interpolations.LinearInterpolation(temp_t,temp)
no3_arr=Interpolations.LinearInterpolation(no3_t,no3)
irr_arr=Interpolations.LinearInterpolation(irr_t,irr)

a_0 = 30;n_0 = 0.009;c_0 = 0.3

solution, results = Kelp.solvekelp(t_i, nd, u_arr, temp_arr, irr_arr, no3_arr, lat, a_0, n_0, c_0)

temp_disp=[];irr_disp=[];n_disp=[];

t_disp=solution.t

pyplot()
plot(layout=grid(2, 3))

plot!(temp_t,temp,sp=1,ylabel="Temperature/degC",legend=true)
plot!(irr_t,irr,sp=2,ylabel="Irradiance/micro mol photons / m^2 / s",legend=true)
plot!(no3_t,no3,sp=3,ylabel="Nitrate/micro mol / L",legend=true)

#sjotun 1993
c_t=[10.016,100.156,130.973,164.872,192.607,224.965,266.568,291.222,343.611,409.097].+7*30
c=[0.312,0.349,0.292,0.195,0.224,0.266,0.281,0.307,0.363,0.373]
n_t=[10.255,96.239,135.681,160.924,195.633,222.454,269.785,288.717,345.514,406.255].+7*30
n=[0.007,0.014,0.021,0.024,0.025,0.022,0.023,0.017,0.012,0.013]
a_t=[16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,342,343,344,345,346,347,348,349,350,351,352,353,354,355,356,357,358,359,360,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375,376,377,378,379,380,381,382,383,384,385,386,387,388,389,390,391,392,393,394,395,396,397,398].+7*30
a=[29.264,29.273,29.273,29.273,29.273,29.273,29.273,29.273,29.273,29.273,29.273,29.273,29.273,29.273,29.273,29.273,29.273,29.273,29.273,29.273,29.273,29.273,29.277,29.366,29.455,29.544,29.633,29.722,29.811,29.9,29.989,30.078,30.167,30.256,30.345,30.434,30.523,30.612,30.701,30.79,30.879,30.968,31.057,31.146,31.24,31.336,31.433,31.529,31.625,31.722,31.818,31.915,32.011,32.107,32.204,32.3,32.397,32.493,32.59,32.686,32.782,32.879,32.975,33.072,33.142,33.207,33.271,33.335,33.399,33.464,33.528,33.592,33.657,33.721,33.785,33.849,33.914,33.978,34.042,34.107,34.171,34.235,34.299,34.364,34.466,34.569,34.672,34.775,34.878,34.981,35.084,35.186,35.289,35.392,35.495,35.598,35.701,35.803,35.906,36.009,36.112,36.215,36.318,36.421,36.523,36.626,36.729,36.832,36.948,37.102,37.256,37.41,37.565,37.719,37.873,38.028,38.182,38.336,38.49,38.645,38.799,38.953,39.107,39.262,39.416,39.551,39.68,39.808,39.937,40.065,40.194,40.322,40.451,40.579,40.708,40.868,41.033,41.198,41.364,41.529,41.694,41.86,42.025,42.19,42.355,42.521,42.617,42.553,42.489,42.424,42.36,42.296,42.231,42.167,42.103,42.039,42,42,42,42,42,42,42,42,42,42,42,42,42.025,42.055,42.085,42.114,42.144,42.174,42.203,42.233,42.263,42.292,42.322,42.352,42.381,42.411,42.441,42.47,42.5,42.53,42.559,42.589,42.619,42.627,42.603,42.579,42.554,42.53,42.506,42.482,42.458,42.434,42.41,42.386,42.362,42.337,42.313,42.289,42.265,42.241,42.217,42.193,42.169,42.145,42.121,42.096,42.072,42.048,42.024,42,42.051,42.103,42.154,42.206,42.257,42.309,42.36,42.411,42.463,42.514,42.566,42.617,42.669,42.72,42.771,42.823,42.874,42.926,42.977,43.028,43.08,43.131,43.183,43.234,43.292,43.369,43.446,43.523,43.601,43.678,43.755,43.832,43.909,43.986,44.063,44.14,44.218,44.295,44.372,44.449,44.526,44.473,44.377,44.28,44.184,44.087,43.991,43.895,43.798,43.702,43.605,43.509,43.413,43.316,43.22,43.123,43.027,42.93,42.834,42.738,42.641,42.584,42.529,42.474,42.419,42.364,42.309,42.253,42.198,42.143,42.088,42.033,41.978,41.923,41.868,41.813,41.758,41.702,41.647,41.592,41.537,41.482,41.427,41.372,41.265,41.15,41.034,40.918,40.802,40.687,40.571,40.455,40.34,40.224,40.108,39.993,39.877,39.761,39.645,39.53,39.425,39.339,39.253,39.167,39.082,38.996,38.91,38.825,38.739,38.653,38.567,38.482,38.396,38.31,38.225,38.182,38.182,38.182,38.182,38.182,38.182,38.182,38.182,38.182,38.182,38.182,38.182,38.182,38.149,38.039,37.928,37.818,37.708,37.598,37.488,37.377,37.267,37.157,37.047,36.937,36.851,36.774,36.697,36.62,36.543,36.466,36.388,36.311,36.234,36.157,36.08,36.003,35.926,35.848,35.771,35.694,35.652,35.717,35.781,35.845,35.91,35.974,36.038,36.102,36.167,36.231,36.295]


plot!(t_disp,results.area,sp=4,xlabel="Month", ylabel="Frond Area/dm^2",ylim=(0,140),label="Model")
plot!(a_t,a,label="Paper",sp=4)
plot!(t_disp,results.nitrogen,sp=5, xlabel="Month", ylabel="Nitrogen reserve/gN/g sw",label="Model")
plot!(n_t,n,sp=5,seriestype=:scatter, label="Sjotun, 1993")
plot!(t_disp,results.carbon,sp=6, xlabel="Month", ylabel="Carbon reserve/gC/g sw",label="Model")
plot!(c_t,c,sp=6,seriestype=:scatter, label="Sjotun, 1993")
xticks=[]
for day in t_i:30:t_i+365
    push!(xticks,Dates.format((Date(2015,1,1)+Dates.Day(day)),"U")[1])
end

plot!(xticks=([t_i:30:t_i+365;],xticks),margin=-3mm)
plot!(top_margin=0mm)
display(display(plot!()))