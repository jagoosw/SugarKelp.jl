# Model parameters
A_0 = 6          # Growth rate adjustment parameter
alpha = 3.75e-5 * 24# photosynthetic efficiency (converted to day^-1)
C_min = 0.01     # Minimal carbon reserve
C_struct = 0.2   # Amount of carbon per unit dry weight of structural mass
gamma = 0.5      # Exudation parameter
epsilon = 0.22   # Frond erosion parameter
I_sat = 200      # Irradiance for maximal photosynthesis
J_max = 1.4e-4 * 24# Maximal nitrate uptake (gN/dm^2/h converted to gN/dm^2/day)
K_A = 0.6        # Structural dry weight per unit area
K_DW = 0.0785    # Dry weight to wet weight ratio of structural mass
K_C = 2.1213     # Mass of carbon reserves per gram carbon
K_N = 2.72       # Mass of nitrogen reserves per gram nitrogen
m_1 = 0.1085     # Growth rate adjustment parameter
m_2 = 0.03       # Growth rate adjustment parameter
mu_max = 0.18    # Maximal area specific growth ratio
N_min = 0.01     # Minimal nitrogen reserve
N_max = 0.022    # Maximal nitrogen reserve
N_struct = 0.01  # Amount of nitrogen per unit dry weight of structural mass
P_1 = 1.22e-3 * 24    # Maximal photosynthetic rate at T = T?P1K converted to day^-1
P_2 = 1.44e-3 * 24    # Maximal photosynthetic rate at T = T?P2K converted to day^-1
a_1 = 0.85       # Photoperiod parameter
a_2 = 0.3        # Photoperiod parameter
R_1 = 2.785e-4 * 24   # Respiration rate at T = TR1, converted to days^-1
R_2 = 5.429e-4 * 24   # Respiration rate at T = TR2, converted to days^-1
T_R1 = 285       # Reference temperature for respiration (K)
T_R2 = 290       # Reference temperature for respiration (K)
T_P1 = 285       # Reference temperature for photosynthesis (K)
T_P2 = 288       # Reference temperature for photosynthesis (K)
T_AP = 1694.4    # Arrhenius temperature for photosynthesis (K)
T_PL = 271
T_PH = 296
T_APH = 25924    # Arrhenius temperature for photosynthesis at high end of range (K)
T_APL = 27774    # Arrhenius temperature for photosynthesis at low end of range (K)
T_AR = 11033     # Arrhenius temperature for respiration (K)
U_0p65 = 0.03    # Current speed at which J = 0.65Jmax (m/s)
K_X = 4 / 1e3          # Nitrate uptake half saturation constant, converted to mmol/L
