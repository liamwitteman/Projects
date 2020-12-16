from co2e_model_inputs import *
import numpy as np
from math import exp


# Calculate double layer initial conditions:
phi_dl_an_0 = phi_an_0 - (phi_elyte_ca_0 - dPhi_j)
phi_dl_ca_0 = phi_elyte_ca_0 - phi_ca_0

P_w = 10**(A-(B/(T+C)))*bar_2_atm
P_ca_0 = P*atm_2_pa
X_k_ca_0 = np.array([1-P_w-1e-12, P_w, 1e-12])
# Integration time is 1000 seconds to achieve steady state.
t_final = 100

"============ INITIALIZE SOLUTION VECTOR ============"
C_k_ca_cl_0 = P_ca_0*X_k_ca_0/R/T
# Initial solution vector:

SV_0 = np.hstack((phi_dl_an_0, phi_dl_ca_0, act_w_aem, C_k_ca_cl_0, C_k_ca_cl_0, 10))

# Create class to point to the correct variable locations in the SV:
class ptr:
    phi_an = 0
    phi_ca = 1

    act_w_aem = 2
    
    # C_k in anode GDL: starts just after phi_dl, is same size as X_k_an:
    C_k_ca_gdl = np.array([3, 4, 5])
    
    # C_k in anode CL: starts just after GDL, is same size as X_k_an:
    C_k_ca_cl = np.array([6, 7, 8])

    # Concentration of water in bipolar membrane junction
    lambda_j = 9
    
# Load inputs and other parameters into 'pars' class:

class pars:
    # Component thicknesses:
    A_cell = A_cell
    dy_gdl_an = dy_gdl_an
    dy_gdl_ca = dy_gdl_ca
    dy_cl_an = dy_cl_an
    dy_cl_ca = dy_cl_ca
    dy_cem = dy_cem
    dy_aem = dy_aem
    dy_j = dy_j

    eps_gdl_an = eps_gdl_an
    eps_gdl_ca = eps_gdl_ca
    eps_cl_an = eps_cl_an
    eps_cl_ca = eps_cl_ca

    n_Brugg_gdl = alpha_Brugg_GDL
    n_Brugg_cl = alpha_Brugg_CL

    d_solid_gdl = d_part_gdl
    d_solid_cl = d_part_cl

    eps_g_dy_Inv_cl = 1/dy_cl_ca/eps_cl_ca

    # Membrane properties
    iec_aem = iec_aem
    iec_cem = iec_cem
    lambda_aem = lambda_aem

    rho_aem = rho_aem*1000
    rho_cem = rho_cem*1000
    Mw_h2o = Mw_h2o/1000

    # Ag surface area per unit geometric area:
    A_fac_Ag = 3*(1-eps_cl_ca)/(d_solid_cl/2)*dy_cl_an*Ag_surf_frac
    # Geometric area per unit double layer area:
    A_fac_dl = Ag_surf_frac/A_fac_Ag
    f_Ag = Ag_surf_frac

    X_k_GDL = X_k_ca_0

    # Equilibrium double layer potentials (V)
    # Assume fixed (for now!)
    dPhi_eq_an = dPhi_eq_an
    dPhi_eq_ca = dPhi_eq_ca
    dPhi_j = dPhi_j

    C_dl_an = C_dl_an
    C_dl_ca = C_dl_ca


    # Butler-Volmer parameters:
    A_k_an = A_k_an
    n_an = n_an
    beta_an = beta_an

    A_k_ca = A_k_ca
    n_ca = n_ca
    beta_ca = beta_ca
    
    # Exchange current density:
    i_o_an = A_k_an*exp(-Ea_an/R/T)
    i_o_ca = A_k_ca*exp(-Ea_ca/R/T)
    # External current density:
    i_array = i_array

    # Electronic and ionic resistivities (ohm-m):
    R_el_gdl_an = 1/sigma_el_nickel/eps_gdl_an
    R_el_cl_an = 1/sigma_el_nickel/eps_cl_an
    R_el_gdl_ca = 1/sigma_el_graphite/eps_gdl_ca
    R_el_cl_ca = 1/sigma_el_graphite/eps_cl_ca

    # Potential drops due to electrical resistance (V)
    dPhi_el_gdl_an = R_el_gdl_an*dy_gdl_an*i_array
    dPhi_el_cl_an = R_el_cl_an*dy_cl_an*i_array
    dPhi_el_gdl_ca = R_el_gdl_ca*dy_gdl_ca*i_array
    dPhi_el_cl_ca = R_el_cl_ca*dy_cl_ca*i_array

    # Operating temperature
    T = T

    # Binary diffusion coefficients
    D_k_g_ca = np.array([(D_co2_co+D_co2_h2o)/2, D_co2_h2o, D_co2_co])

    mu_g_ca = mu_g_ca
    nu_k_ca = nu_k_ca

    # Saturation pressure of water from Antoine's equation
    P_w = P_w

    bar_2_atm = bar_2_atm
    atm_2_pa = atm_2_pa

