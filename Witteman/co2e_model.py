""" CO2E
 model."""
from scipy.integrate import solve_ivp
import numpy as np
from co2e_function import residual
from co2e_model_init import SV_0, t_final, pars, ptr
from matplotlib import pyplot as plt
from math import exp

# Constants
F = 96485
R = 8.3145

i_array = pars.i_array
V_cell = np.zeros_like(i_array)

lambda_j = np.zeros_like(i_array)

act_w_aem = np.zeros_like(i_array)
act_avg_aem = np.zeros_like(i_array)
cond_io_aem = np.zeros_like(i_array)
dPhi_io_aem = np.zeros_like(i_array)

act_w_cem = np.zeros_like(i_array)
lambda_cem = np.zeros_like(i_array)
lambda_avg_cem = np.zeros_like(i_array)
cond_io_cem = np.zeros_like(i_array)
dPhi_io_cem = np.zeros_like(i_array)


time_span = np.array([0,t_final])
SV = np.zeros((SV_0.size,i_array.size))
for j,i_ext in enumerate(i_array):
    solution = solve_ivp(lambda t, y: residual(t, y, pars, ptr, i_ext), 
        time_span, SV_0, rtol=1e-9, atol=1e-8, method='BDF')
    
    SV[:,j] = solution.y[:,-1]
    SV_0 = SV[:,j]
    lambda_j[j] = SV[ptr.lambda_j,j]
    # Post processing AEM
    
    act_avg_aem[j] = lambda_j[j]/pars.lambda_aem
    cond_io_aem[j] = 100*((0.0018*act_avg_aem[j]**4-0.0036*act_avg_aem[j]**3+0.0025*act_avg_aem[j]**2-0.00059*act_avg_aem[j]+3.71e-5)*(pars.T-333) + (0.011*act_avg_aem[j]**4-0.0064*act_avg_aem[j]**3+0.022*act_avg_aem[j]**2+0.00044*act_avg_aem[j]+8.19e-5))
    dPhi_io_aem[j] = (1/cond_io_aem[j])*pars.dy_aem*i_ext   

    # Post processing CEM
    act_w_cem[j] = SV[7,j]*R*pars.T/pars.atm_2_pa/pars.P_w
    if act_w_cem[j] <= 1:
        lambda_cem[j] = 0.043+17.18*act_w_cem[j]-39.85*act_w_cem[j]**2+36.0*act_w_cem[j]**3
    else:
        lambda_cem[j] = 14+1.4*(act_w_cem[j]-1)
    lambda_avg_cem[j] = (lambda_j[j]+lambda_cem[j])/2
    cond_io_cem[j] = 100*((0.005193*lambda_avg_cem[j])*exp(1268*(1/303.15-1/pars.T)))
    dPhi_io_cem[j] = (1/cond_io_cem[j])*pars.dy_cem*i_ext


V_an = SV[ptr.phi_an,:] + SV[ptr.phi_ca,:] - pars.dPhi_j
V_el = pars.dPhi_el_gdl_an + pars.dPhi_el_cl_an + pars.dPhi_el_gdl_ca + pars.dPhi_el_cl_ca
V_io = dPhi_io_aem + dPhi_io_cem

V_cell = V_an + V_el + V_io
i_ext = i_array/1e4

i_file = open("i_data.txt","w")
np.savetxt(i_file, i_ext, delimiter=",")

v_file = open("v_data.txt","w")
np.savetxt(v_file, SV[6,:]*R*pars.T/pars.atm_2_pa, delimiter=",")

fig, ax = plt.subplots()
plt.plot(i_ext,V_cell) 

plt.xlabel('Current Density (A/cm$^2$)',fontsize=14)
plt.ylabel('Applied Potential (V)',fontsize=14)
plt.legend(['Model','Experimental'],frameon=False)
plt.show()


for var in SV[ptr.C_k_ca_cl,:]:
    plt.plot(i_ext,var*R*pars.T/pars.atm_2_pa)
    
plt.legend(['CO2','H2O','CO'])
plt.xlabel('Current Density (A/cm$^2$)',fontsize=14)
plt.ylabel('Species partial pressures in CL (atm)',fontsize=14)
plt.show()

plt.plot(i_ext,SV[ptr.lambda_j,:])
    
plt.xlabel('Current Density (A/cm$^2$)',fontsize=14)
plt.ylabel('Water hydration in BPM junction (mol H$_2$O/mol SO$_4$-)',fontsize=14)
plt.show()
