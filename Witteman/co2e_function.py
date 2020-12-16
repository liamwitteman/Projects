import numpy as np
from math import exp

# Constants
F = 96485
R = 8.3145

def residual(t, SV, pars, ptr, i_ext):
    dSV_dt = np.zeros_like(SV)

    # Calculations to get activity of each species
    # GDL gas phase:
    C_k_ca_gdl = SV[ptr.C_k_ca_gdl]
    # Catalyst layer gas phase
    C_k_ca_cl = SV[ptr.C_k_ca_cl]

    P_k_ca_cl = C_k_ca_cl*R*pars.T
    P_tot = sum(P_k_ca_cl)

    act_co2 = P_k_ca_cl[0]/P_tot
    act_h2o = P_k_ca_cl[1]/P_tot
    act_co = P_k_ca_cl[2]/P_tot
    act_w_cem = P_k_ca_cl[1]/pars.atm_2_pa/pars.P_w
    act_w_aem = SV[ptr.act_w_aem]
    

    # Anode Interface:
    eta_an = SV[ptr.phi_an] - pars.dPhi_eq_an

    i_Far_an = pars.i_o_an*(exp(pars.n_an*F*pars.beta_an*eta_an/R/pars.T)
                            - exp(-pars.n_an*F*(1-pars.beta_an)*eta_an/R/pars.T))
    i_dl_an = i_ext*pars.A_fac_dl - i_Far_an*pars.f_Ag
    dSV_dt[ptr.phi_an] = i_dl_an/pars.C_dl_an

    # Cathode Interface:
    eta_ca = SV[ptr.phi_ca] - (pars.dPhi_eq_ca - R*pars.T/pars.n_an/F*np.log(act_h2o*act_co/act_co2))
    i_Far_ca = pars.i_o_ca*(exp(pars.n_ca*F*pars.beta_ca*eta_ca/R/pars.T)
                            - exp(-pars.n_ca*F*(1-pars.beta_ca)*eta_ca/R/pars.T))
    i_dl_ca = i_ext*pars.A_fac_dl - i_Far_ca*pars.f_Ag
    dSV_dt[ptr.phi_ca] = i_dl_ca/pars.C_dl_ca

    dSV_dt[ptr.act_w_aem] = 0


   
    s1 = {'C_k': C_k_ca_gdl, 'dy':pars.dy_gdl_ca, 'eps_g':pars.eps_gdl_ca, 
        'n_Brugg':pars.n_Brugg_gdl, 'd_solid':pars.d_solid_gdl}
    s2 = {'C_k': C_k_ca_cl, 'dy':pars.dy_cl_ca, 'eps_g':pars.eps_cl_ca,
        'n_Brugg':pars.n_Brugg_cl, 'd_solid':pars.d_solid_cl}
    gas_props = {'T':pars.T, 'D_k':pars.D_k_g_ca, 'mu':pars.mu_g_ca}
    
    N_k_i = co2e_gas_flux(s1, s2, gas_props)

    # Molar production rates due to Faradaic current:
    sdot_k = i_Far_ca*pars.nu_k_ca/pars.n_ca/F
    
    # Change in catalyst layer gas phase mole fractions:
    dCk_dt = (N_k_i + sdot_k)*pars.eps_g_dy_Inv_cl

    # water flux calculations through cem
    if act_w_cem <= 1:
        lambda_cem = 0.043+17.18*act_w_cem-39.85*act_w_cem**2+36.0*act_w_cem**3
    else:
        lambda_cem = 14+1.4*(act_w_cem-1)
    
    lambda_avg_cem = (lambda_cem+SV[ptr.lambda_j])/2
    D_w_cem = (2.563-0.33*lambda_avg_cem+0.0264*lambda_avg_cem**2-0.000671*lambda_avg_cem**3)*1e-10*exp(2417*(1/303.15-1/pars.T))
    n_d_cem = 2.5*lambda_avg_cem/22
    f_w_cl_ca = lambda_cem*pars.Mw_h2o/pars.iec_cem
    f_w_int_ca = SV[ptr.lambda_j]*pars.Mw_h2o/pars.iec_cem
    f_w_avg_cem = (f_w_cl_ca+f_w_int_ca)/2
    J_w_cem = -pars.rho_cem*D_w_cem/pars.Mw_h2o/(1+f_w_avg_cem)*(f_w_cl_ca-f_w_int_ca)/pars.dy_cem+n_d_cem*i_ext/F
   
    N_cem_i = np.array([0, J_w_cem, 0])

    dCk_dt = dCk_dt + N_cem_i/pars.dy_cem
    dSV_dt[ptr.C_k_ca_cl] = dCk_dt

    # water flux calculations through cem
    lambda_avg_aem = (pars.lambda_aem+SV[ptr.lambda_j])/2
    D_w_aem = 6e-10
    n_d_aem = 1
    f_w_cl_an = pars.lambda_aem*pars.Mw_h2o/pars.iec_aem
    f_w_int_an = SV[ptr.lambda_j]*pars.Mw_h2o/pars.iec_aem
    f_w_avg_aem = (f_w_cl_an+f_w_int_an)/2
    J_w_aem = -pars.rho_aem*D_w_aem/pars.Mw_h2o/(1+f_w_avg_aem)*(f_w_cl_an-f_w_int_an)/pars.dy_aem+n_d_aem*i_ext/F

    dSV_dt[ptr.lambda_j] = (-J_w_cem - J_w_aem - i_ext/F)/1000/pars.dy_j

    return dSV_dt


def co2e_gas_flux(node1, node2, gas_props):
    N_k  = np.zeros_like(node1['C_k'])

    f1 = node1['dy']/(node1['dy'] + node2['dy'])
    f2 = 1-f1

    C_int = f1*node1['C_k'] + f2*node2['C_k']

    X_k_1 = node1['C_k']/np.sum(node1['C_k'])
    X_k_2 = node2['C_k']/np.sum(node2['C_k'])
    X_k_int = f1*X_k_1 + f2*X_k_2


    P_1 = np.sum(node1['C_k'])*R*gas_props['T']
    P_2 = np.sum(node2['C_k'])*R*gas_props['T']
    # print(P_1, P_2)

    eps_g = f1*node1['eps_g'] + f2*node2['eps_g']
    tau_fac = (f1*node1['eps_g']**node1['n_Brugg'] 
        + f2*node2['eps_g']**node2['n_Brugg'])
    D_k_eff = eps_g*gas_props['D_k']/tau_fac
    # print(f2, f1)
    # print(tau_fac)
    # print(f1*node1['eps_g']**node1['n_Brugg'])
    
    d_part = f1*node1['d_solid'] + f2*node2['d_solid']
    K_g = eps_g**3*d_part**2*tau_fac**(-2)*(1-eps_g)**(-2)/72

    dY = 0.5*(node1['dy'] + node2['dy'])

    V_conv = -K_g*(P_2 - P_1)/dY/gas_props['mu']
    V_k_diff = -D_k_eff*(X_k_2 - X_k_1)/dY/X_k_int
    V_k  = V_conv + V_k_diff
    N_k = C_int*X_k_int*V_k

    return N_k