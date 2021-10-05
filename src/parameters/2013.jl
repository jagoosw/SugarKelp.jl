const A_0 = 4.5#6          # Growth rate adjustment parameter
const α = 4.15e-5 * 24 * 10^6 / (24 * 60 * 60)#3.75e-5 * 24 * 10^6 / (24 * 60 * 60) # photosynthetic efficiency 
const C_min = 0.01     # Minimal carbon reserve
const C_struct = 0.2   # Amount of carbon per unit dry weight of structural mass
const γ = 0.5# Exudation parameter
const ϵ = 0.22   # Frond erosion parameter
const I_sat = 90 * 24*60*60/(10^6)#200 * 24 * 60 * 60 / (10^6)  # Irradiance for maximal photosynthesis
const J_max = 1.4e-4 * 24# Maximal nitrate uptake (gN/dm^2/h converted to gN/dm^2/day)
const K_A = .5#0.6        # Structural dry weight per unit area
const K_DW = 0.0785    # Dry weight to wet weight ratio of structural mass
const K_C = 2.1213     # Mass of carbon reserves per gram carbon
const K_N = 2.72       # Mass of nitrogen reserves per gram nitrogen
const N_min = 0.0126#0.01     # Minimal nitrogen reserve
const N_max = 0.0216#0.022    # Maximal nitrogen reserve
const m_2 = 0.039/(2*(1-N_min/N_max))#0.036# 0.03       # Growth rate adjustment parameter
const m_1 = 0.18/(2*(1-N_min/N_max))-m_2# 0.1085     # Growth rate adjustment parameter
const μ_max = 0.18    # Maximal area specific growth ratio
const N_struct = 0.0146#0.01  # Amount of nitrogen per unit dry weight of structural mass
const P_1 = 1.22e-3 * 24    # Maximal photosynthetic rate at T = T?P1K converted to day^-1
const P_2 = 1.3e-3 * 24    # Maximal photosynthetic rate at T = T?P2K converted to day^-1
const a_1 = 0.85       # Photoperiod parameter
const a_2 = 0.3        # Photoperiod parameter
const R_1 =  2.785e-4 * 24   # Respiration rate at T = TR1,
const R_2 = 5.429e-4 * 24   # Respiration rate at T = TR2,
const T_R1 = 285       # Reference temperature for respiration (K)
const T_R2 = 290       # Reference temperature for respiration (K)
const T_P1 = 285       # Reference temperature for photosynthesis (K)
const T_P2 = 288       # Reference temperature for photosynthesis (K)
const T_AP = (1 / T_P1 - 1 / T_P2)^(-1) * log(P_2 / P_1)# 1694.4    # Arrhenius temperature for photosynthesis (K)
const T_PL = 271
const T_PH = 296
const T_APH = 1414.87# 1516.7#25924    # Arrhenius temperature for photosynthesis at high end of range (K)
const T_APL = 4547.89# 4388.5#27774    # Arrhenius temperature for photosynthesis at low end of range (K)
const T_AR = (1 / T_R1 - 1 / T_R2)^(-1) * log(R_2 / R_1)     # Arrhenius temperature for respiration (K)
const U_0p65 = 0.03    # Current speed at which J = 0.65Jmax (m/s)
const K_X = 4          # Nitrate uptake half saturation constant, converted to mmol/m^3

const R_A=1.11e-4*24
const R_B=5.57e-5*24