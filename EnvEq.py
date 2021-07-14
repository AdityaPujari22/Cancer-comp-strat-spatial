# All equations for modelling tumour environment

import scipy
from scipy.integrate import solve_ivp,odeint,ode
import seaborn as sns
import pandas as pd

# Differential equation maker
def enveq(t, x , p , mu , lam , r , K , delta , rho , K_m, lim , D, M):
        
    # Seperating x into lists of variables, of length M
    o2 = []
    test = []
    Tpro = []
    Tpos = []
    Tneg = []
    for i in range(len(x)):
        if i % 5 == 0:
            Tpos.append(x[i])
                
        if i % 5 == 1:
            Tpro.append(x[i])
            
        if i % 5 == 2:
            Tneg.append(x[i])
                
        if i % 5 == 3:
            o2.append(x[i])
                
        if i % 5 == 4:
            test.append(x[i])
     
    do2 = []
    dtest = []
    dTpro = []
    dTpos = []
    dTneg = []
        
    for i in range(M):
        
        if i == 0: #If first compartment, resource only diffuses to second compartment
            #Equation for oxygen: constant production, uptake by all 3 cells, decay, diffusion
            d3 = p[0][i] - mu[0,0] * Tpos[i] - mu[0,1] * Tpro[i] - mu[0,2] * Tneg[i] - lam[0] * o2[i] + D[0] * (o2[i + 1] - 2 * o2[i])
            do2.append(d3)
            
            #Equation for testosterone: production by Tp, uptake by all Tp, T+, decay, diffusion
            d4 = p[1] * Tpro[i] - mu[1,0] * Tpos[i] - mu[1,1] * Tpro[i] - lam[1] * test[i] + D[1] * (test[i + 1] - 2 * test[i])
            dtest.append(d4)
        
        if i == M-1: #If last compartment, resource only diffuses to second-to-last compartment
            #Equation for oxygen: constant production, uptake by all 3 cells, decay, diffusion
            d3 = p[0][i] - mu[0,0] * Tpos[i] - mu[0,1] * Tpro[i] - mu[0,2] * Tneg[i] - lam[0] * o2[i] + D[0] * (o2[i - 1] - 2 * o2[i]) 
            do2.append(d3)
            
            #Equation for testosterone: production by Tp, uptake by all Tp, T+, decay, diffusion
            d4 = p[1] * Tpro[i] - mu[1,0] * Tpos[i] - mu[1,1] * Tpro[i] - lam[1] * test[i] + D[1] * (test[i - 1] - 2 * test[i])
            dtest.append(d4)
            
        elif 0 < i < M-1: #Otherwise, resource diffuses in either direction
            #Equation for oxygen: constant production, uptake by all 3 cells, decay, diffusion
            d3 = p[0][i] - mu[0,0] * Tpos[i] - mu[0,1] * Tpro[i] - mu[0,2] * Tneg[i] - lam[0] * o2[i] + D[0] * (o2[i - 1] + o2[i + 1] - 2 * o2[i])
            do2.append(d3)
            
            #Equation for testosterone: production by Tp, uptake by all Tp, T+, decay, diffusion
            d4 = p[1] * Tpro[i] - mu[1,0] * Tpos[i] - mu[1,1] * Tpro[i] - lam[1] * test[i] + D[1] * (test[i - 1] + test[i + 1] - 2 * test[i])
            dtest.append(d4)
    
        #Equation for T+
        d0 = r[0] * Tpos[i] * (1 - ((Tpos[i] + Tpro[i] + Tneg[i])/(K + rho[0] * f_res(o2[i] , lim[0,0])*f_res(test[i],lim[1,0])))) - delta[0] * Tpos[i] 
        
        if Tpos[i] > K_m[0]: #If compartment exceeds maximum capacity for cells
            ovr = Tpos[i] % K_m[0]
            Tpos[i] = K_m[0]
            if i == 0: # If first compartment, sends cells to the second compartment
                Tpos[i + 1] += ovr  
                
            if i == M-1: #If last compartment, sends cells to M-1th compartment
                Tpos[i - 1] += ovr
                
            if 0 < i < M-1: #Otherwise equally divides overflow and sends to compartments before and after
                Tpos[i + 1] = 0.5*ovr
                Tpos[i - 1] = 0.5*ovr
           
        #Equation for Tp 
        d1 = r[1] * Tpro[i] * (1 - ((Tpos[i] + Tpro[i] + Tneg[i])/(K + rho[1] * f_res(o2[i], lim[0,1]) * f_res(test[i],lim[1,1])))) - delta[1] * Tpro[i]
        
        if Tpro[i] > K_m[1]: #If compartment exceeds maximum capacity for cells
            ovr = Tpro[i] % K_m[1]
            Tpro[i] = K_m[1]
            if i == 0: # If first compartment, sends cells to the second compartment
                Tpro[i + 1] += ovr  
                
            if i == M-1: #If last compartment, sends cells to M-1th compartment
                Tpro[i - 1] += ovr
               
            if 0 < i < M-1: #Otherwise equally divides overflow and sends to compartments before and after
                Tpro[i + 1] = 0.5*ovr
                Tpro[i - 1] = 0.5*ovr
                       
        #Equation for T-    
        d2 = r[2] * Tneg[i] * (1 - ((Tpos[i] + Tpro[i] + Tneg[i])/(K + rho[2] * f_res(o2[i],lim[0,2])))) - delta[2] * Tneg[i]
        
        if Tneg[i] > K_m[2]: #If compartment exceeds maximum capacity for cells
            ovr = Tneg[i] % K_m[2]
            Tneg[i] = K_m[2]
            if i == 0: # If first compartment, sends cells to the second compartment
                Tneg[i + 1] += ovr  
               
            if i == M-1: #If last compartment, sends cells to M-1th compartment
                Tneg[i - 1] += ovr
                
            if 0 < i < M-1: #Otherwise equally divides overflow and sends to compartments before and after
                Tneg[i + 1] = 0.5*ovr
                Tneg[i - 1] = 0.5*ovr                
        
        dTpro.append(d1)
        dTpos.append(d0)
        dTneg.append(d2)
     
    # Returns the list with dx_i_m/dt at that time
    dx = []
    for i in range(M):
        dx.append(dTpos[i])
        dx.append(dTpro[i])
        dx.append(dTneg[i])
        dx.append(do2[i])
        dx.append(dtest[i])

    return dx

#Resource Uptake function
def f_res(res, k):
    ll = k[0]
    ul = k[1]
    
    if res <= ll:
        return 0
    if res >= ul:
        return 1
    else:
        return (res - ll)/(ul - ll)    
    
#Function for solving system
def solve_eq(t_max, dt, y0, p, mu, lam, r, K, delta, rho, K_m, lim, D, M):
    
    #Timeseries arrays
    t0 = 0
    t = [t0]
    y_0 = np.reshape(y0, (M, 5))
    y = [y_0]
    y_t = y0
    
    f_ode = ode(lambda t,y: enveq(t , y , p , mu , lam , r , K , delta , rho , K_m, lim , D, M)).set_integrator('lsoda')
    f_ode.set_initial_value(y0, t0)    
    while f_ode.t < t_max:
        t.append(f_ode.t+dt)
        y_t = f_ode.integrate(f_ode.t+dt)
                      
        for i in range(len(y_t)):
            if y_t[i] < 0: # If any value goes negative
                y_t[i] = 0 # Sets neg values to 0
                f_ode.set_initial_value(y_t,f_ode.t) #if anomaly detected, reset initial conditions to zero set values at that time t

        y_t = np.reshape(y_t, (M, 5))
        
        for i in range(len(y_t)):
            if np.logical_and(y_t[i][0:3]>0 , y_t[i][0:3]<1).any(): # If cell numbers drop below 1
                y_t[i][0:3] = np.where(y_t[i][0:3] >= 1, y_t[i][0:3], np.zeros(np.shape(y_t[i][0:3]))) #Set cell number to 0                
        y.append(y_t)        
        if not f_ode.successful(): #Check if integration fails
            print('Failure @',f_ode.t)
                
    result= [t, y]    
    return result