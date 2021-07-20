#All input parameters 

import numpy as np


# Number of compartments ( > 1)
M = 3

#Time parameters of simulation
t_max = 1000*24*60 #maximum (mins)
dt = 60  #step (mins)
t_steps = t_max / dt

#Initial Values

##cells
y0_Tpos = np.zeros(M) 
y0_Tpro = np.zeros(M) 
y0_Tneg = np.zeros(M) 

##resource
y0_o2 = np.ones(M) * 0.5 #(prop)
y0_test= np.zeros(M) #(prop)


#Production rates of resources

p_o2 = np.ones(M) * 0.11 #(prop/min)
p_test = 5E-7 #by Tp cells (prop/min/cell)


#Uptake rate of resources by cells

##oxygen
mu_o2Tpos = 1.63E-6 #(prop/min/cell)
mu_o2Tpro = 1.63E-6 #(prop/min/cell)
mu_o2Tneg = 1.04E-6 #(prop/min/cell)

##testosterone
mu_testTpos = 2.34E-8 #(prop/min/cell)
mu_testTpro = 6E-8 #(prop/min/cell)

#Decay rates of resources
lam_o2 = 0.1 #(1/min)
lam_test = 0.004 #(1/min)

#Doubling time
t_DTpos = 34*60 #(min)
t_DTpro = 40*60 #(min)
t_DTneg = 25*60 #(min)

#Growth rate
r_Tpos = 2.84E-3 #(/min)
r_Tpro = 2.79E-3 #(/min)
r_Tneg = 6.23E-4 #(/min)

#Death rate
delta_Tpos = 2.5E-3 #(/min)
delta_Tpro = 2.5E-3 #(/min)
delta_Tneg = 1.6E-4 #(/min)

#Min Carrying capacity
K = 1

#Environmental Carrying capacity (Scaling Factor)
rho_Tpos = 8.35E4
rho_Tpro = 9.62E4
rho_Tneg = 1.34E4

# Maximum Capacity per compartment
K_Tpos = (K + rho_Tpos)/M
K_Tpro = (K + rho_Tpro)/M
K_Tneg = (K + rho_Tneg)/M

#Resource limits

##Oxygen
###T+
l_lim_o2Tpos = 0 #lower-threshold(prop)
u_lim_o2Tpos = 1.1 #upper-saturation(prop)

###Tp
l_lim_o2Tpro = 0 #(prop)
u_lim_o2Tpro = 1.1 #(prop)

###T-
l_lim_o2Tneg = 0 #(prop)
u_lim_o2Tneg = 1.1 #(prop)

##Testosterone
###T+
l_lim_testTpos = 0 #(prop)
u_lim_testTpos = 0.1 #(prop)

###Tp
l_lim_testTpro = 0 #(prop)
u_lim_testTpro = 0.1 #(prop)

###Diffusion Coefficients

# Sources and values: 
# https://www.nature.com/articles/2151173a0 - 3.6E-5 in rat liver
# https://pubmed.ncbi.nlm.nih.gov/563582/ - 1.75E-5 in rat carcinoma
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1204697/ - 1.64E-5 in rat liver

D_o2 = 60 * 1.7E-5 # (cm^2/min)

# https://www.sciencedirect.com/science/article/abs/pii/S0022354915503316 - 1.9E-6 in pig intestinal mucous
D_test = 60 * 1.9E-6 # (cm^2/min)



# Creating initial values dataset, by parsing input into numpy arrays

y0 = []
for i in range(M):
    y0.append(y0_Tpos[i])
    y0.append(y0_Tpro[i])
    y0.append(y0_Tneg[i])
    y0.append(y0_o2[i])
    y0.append(y0_test[i])

p = [p_o2, p_test]

mu = np.array([[mu_o2Tpos,mu_o2Tpro,mu_o2Tneg],[mu_testTpos,mu_testTpro,0]])

lam = np.array([lam_o2,lam_test])

t_D = np.array([t_DTpos,t_DTpro,t_DTneg])

r = np.array([r_Tpos,r_Tpro,r_Tneg])

delta = np.array([delta_Tpos,delta_Tpro,delta_Tneg])

rho = np.array([rho_Tpos,rho_Tpro,rho_Tneg])

K_m = np.array([K_Tpos, K_Tpro, K_Tneg])

lim = np.array([[[l_lim_o2Tpos,u_lim_o2Tpos],[l_lim_o2Tpro,u_lim_o2Tpro],[l_lim_o2Tneg,u_lim_o2Tneg]],[[l_lim_testTpos,u_lim_testTpos],[l_lim_testTpro,u_lim_testTpro],[0,0]]],dtype=np.float64)

D = np.array([D_o2, D_test])

dataset = [t_max, dt, y0, p, mu, lam, r, K, delta, rho, K_m, lim, D, M]