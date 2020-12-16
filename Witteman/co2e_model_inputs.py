import numpy as np
from math import exp
from co2e_model_diff import Dab

# Inputs:
i_array = np.concatenate((np.linspace(0, 500, 50), np.linspace(510, 11000, 100)), axis=None)
T = 25+273.15 # Starting operating temnperature, K
R = 8.3145 # Ideal gas constant, J/mol-K
P = 1 # Starting operating pressure, atm

" Microstructure and geometry "
dy_gdl_an = 1600e-6 # Anode gas diffusion layer thickness (of nickel foam), m
dy_gdl_ca = 365e-6 # Cathode gas diffusion layer thickness, m
dy_cl_an = 5e-6  # Anode catalyst layer thickness, m
dy_cl_ca = 5e-6 # Cathode catalyst layer thickness, m
dy_cem = 25e-6 # Cation exchange membrane thickness, m
dy_aem = 25e-6 # Anion exchange membrane thickness, m
dy_j = 4e-9 # Bipolar membrane junction thickness, m

eps_gdl_an = 0.95 # Porosity of anode diffusion media (nickel foam)
eps_gdl_ca = 0.78 # Porosity of cathode gas diffusion layer (carbon cloth)
eps_cl_an = 0.95 # Porosity of anode catalyst layer
eps_cl_ca = 0.5 # Porosity of cathode catalyst layer

# Exponent in the Bruggeman correlatin: tau_fac = eps^alpha
alpha_Brugg_GDL = -1
alpha_Brugg_CL = -0.5

# Carbon particle diameter:
d_part_gdl = 2*7.33e-7 # m
d_part_cl = 2*5e-8 # m
Ag_surf_frac = 0.1 # Fraction of carbon surface covered by Ag in CL.
A_cell = 4/100/100 # Active area of cell, m2

" Electrochemical parameters "
# Initial voltages (used to calculate phi_dl guesses)
phi_an_0 = 2.0 # V
phi_elyte_an_0 = 1.43  # V
phi_elyte_ca_0 = 0.6 # V
phi_ca_0 = 0 # V
dPhi_j = -0.828 # V

# Equilibrium double layer potential (phi_electrode - phi_elyte)
dPhi_eq_an = -0.44
dPhi_eq_ca = -0.11

C_dl_an = 1e2 # Double layer capacitance, F/m2 of dl interface
C_dl_ca = 1e2 # Double layer capacitance, F/m2 of dl interface

A_k_an = 1.23e-3  # Anodic exchange current density, A/m2
Ea_an = 25000 # Anodic activation energy, J/mol
n_an = 2 # Charge equivalents transfered to anode
beta_an = 0.5 # Symmetry parameter


A_k_ca = 7.25e8 # Cathodic exchange current density, A/m2
Ea_ca = 100000 # Cathodic activation energy, J/mol
n_ca = 2 # Charge equivalents transfered to cathode
beta_ca = 0.5 # Symmetry parameter

nu_k_ca = np.array([-1., 1., 1.])

" Material properties "
# Anode
density_nickel = 8900 # mass density, kg/m3
# Thermal conductivity from: 
# https://tfaws.nasa.gov/wp-content/uploads/TFAWS18-PT-11.pdf
lambda_cond_an = 1.4 #W/m-K
# conductivity taken from:
# https://www.tibtech.com/conductivite.php?lang=en_US
sigma_el_nickel = 14.3e6 #S/m

# Anion exchange membrane
# Anion water content from Adam Weber paper
act_w_aem = 1 # Activity of water on AEM side (unity since liquid water)
iec_aem = 1.7 # mol/kg
rho_aem = 1.2 # g/ml
lambda_aem = 19

# Cation exchange membrane
act_w_cem = 1 # Activity of water on CEM side (unity since liquid water)
iec_cem = 0.93 # mol/kg
rho_cem = 1.58 # g/ml

# Cathode
density_graphite = 2260 # mass density, kg/m3
# Cp from https://webbook.nist.gov/cgi/cbook.cgi?ID=C7782425&Mask=2
Cp_graphite = 691.67 #J/kg-K 
# conductivity taken as an average from:
# https://en.m.wikipedia.org/wiki/Electrical_resistivity_and_conductivity#Resistivity_of_various_materials
sigma_el_graphite = 2e4 #S/m

# Antoine's constants for water vapor
A = 4.6543
B = 1435.264
C = -64.848

" Transport properties "
Pc_co2 = 72.8 
Pc_h2o = 218.3
Pc_co = 34.5

Tc_co2 = 304.2
Tc_h2o = 647.3
Tc_co = 132.9

Mw_co2 = 44.01 
Mw_h2o = 18.02
Mw_co = 28.01

a_wo_h2o = 2.745e-4
b_wo_h2o = 1.823

a_w_h2o = 3.64e-4
b_w_h2o = 2.334

wo_h2o = {'a': a_wo_h2o, 'b': b_wo_h2o}
w_h2o = {'a': a_w_h2o, 'b': b_w_h2o}

co2 = {'Pc_a': Pc_co2, 'Tc_a': Tc_co2, 'Mw_a': Mw_co2}
h2o = {'Pc_b': Pc_h2o, 'Tc_b': Tc_h2o, 'Mw_b': Mw_h2o}
co = {'Pc_b': Pc_co, 'Tc_b': Tc_co, 'Mw_b': Mw_co}

D_co2_h2o = Dab(co2,h2o,w_h2o, P, T) # binary diffusion coefficient of h2o in co2, m2/s
D_co2_co = Dab(co2,co,wo_h2o, P, T) # binary diffusion coefficient of co in co2, m2/s
mu_g_ca = 16.61e-6 # dynamic viscosity of co2, Pa*s
" Conversions "
bar_2_atm = 0.986923
atm_2_pa = 101325

