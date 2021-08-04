# All equations for modelling tumour environment

from scipy.integrate import solve_ivp,odeint,ode
import pandas as pd
import numpy as np
from pathlib import Path

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
def solve_eq(t_max, dt, y0, p, mu, lam, r, K, delta, rho, K_m, lim, D, M, f_name):
    
    T = np.arange(0, t_max, dt)    
    result = odeint(enveq, y0, T, args = (p, mu, lam, r, K, delta, rho, K_m, lim, D, M))
    result = [np.reshape(i, (M, 5)) for i in result]
    
    t_steps = []
    for i in T:
        for j in range(M):
            t_steps.append(i/(24 * 60))  # Time steps in days
    data = {"Time": t_steps}    
 
     # Returns time-series data as dataframes of values    
    tpos = []
    tpro = []
    tneg = []
    o2 = []
    test = []
    comps = []
    for i in range(len(result)):
        for m in range(M):
            tpos.append(result[i][m, 0])
            tpro.append(result[i][m, 1])
            tneg.append(result[i][m, 2])
            o2.append(result[i][m, 3])
            test.append(result[i][m, 4])
            comps.append(m)
    
    dict0 = {"Tpos" : tpos, "Tpro" : tpro, "Tneg" : tneg, "O2" : o2, "Test" : test, "Compartment" : comps}
    data.update(dict0)
   
    cell_tot = []
    for i in range(len(tpos)):
        cell_tot.append(tpos[i] + tpro[i] + tneg[i]) # Total Cell count per compartment
    dict1 = {"Total Cells" : cell_tot}
    data.update(dict1)
                           
    tot_cell = []
    for i in range(len(result)): # Counting total cell population across all compartments
        k = 0
        for m in range(M):
            k += result[i][m, 0] + result[i][m, 1] + result[i][m, 2]
            
        tot_cell.append(k)
    t_burd = []
    for i in tot_cell:
        for m in range(M):
            t_burd.append(i) # Assigning tumour burden numbers to each time point
                        
    dict3 = {"Tumour Burden" : t_burd}
    data.update(dict3)
    df = pd.DataFrame(data)    
    
    output_file = f_name + ".csv"                     
    output_dir = Path('../../raw_output')
    output_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_dir / output_file)    
    
    return df


# Differential equation maker
def enveq(x, t, p, mu, lam, r, K, delta, rho, K_m, lim, D, M):
        
    # Seperating x into lists of variables, on length M
    o2 = []
    test = []
    Tpro = []
    Tpos = []
    Tneg = []
    for i in range(len(x)):
        if x[i] < 0:
            x[i] = 0
        if i % 5 == 0:
            if 0 < x[i] < 1:
                x[i] = 0
            Tpos.append(x[i])
                
        if i % 5 == 1:
            if 0 < x[i] < 1:
                x[i] = 0
            Tpro.append(x[i])
            
        if i % 5 == 2:
            if 0 < x[i] < 1:
                x[i] = 0
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
        
    for i in range(int(len(x)/5)):
        
        if i == 0:
            #Equation for oxygen: constant production, uptake by all 3 cells, decay
            d3 = p[0][i] - mu[0,0] * Tpos[i] - mu[0,1] * Tpro[i] - mu[0,2] * Tneg[i] - lam[0] * o2[i] + D[0] * (o2[i + 1] - 2 * o2[i])
            do2.append(d3)
            
            #Equation for testosterone: production by Tp, uptake by all Tp, T+, decay
            d4 = p[1] * Tpro[i] - mu[1,0] * Tpos[i] - mu[1,1] * Tpro[i] - lam[1] * test[i] + D[1] * (test[i + 1] - 2 * test[i])
            dtest.append(d4)
        
        if i == M-1:
            #Equation for oxygen: constant production, uptake by all 3 cells, decay
            d3 = p[0][i] - mu[0,0] * Tpos[i] - mu[0,1] * Tpro[i] - mu[0,2] * Tneg[i] - lam[0] * o2[i] + D[0] * (o2[i - 1] - 2 * o2[i]) 
            do2.append(d3)
            
            #Equation for testosterone: production by Tp, uptake by all Tp, T+, decay
            d4 = p[1] * Tpro[i] - mu[1,0] * Tpos[i] - mu[1,1] * Tpro[i] - lam[1] * test[i] + D[1] * (test[i - 1] - 2 * test[i])
            dtest.append(d4)
        elif 0 < i < M-1:
            #Equation for oxygen: constant production, uptake by all 3 cells, decay
            d3 = p[0][i] - mu[0,0] * Tpos[i] - mu[0,1] * Tpro[i] - mu[0,2] * Tneg[i] - lam[0] * o2[i] + D[0] * (o2[i - 1] + o2[i + 1] - 2 * o2[i])
            do2.append(d3)
            
            #Equation for testosterone: production by Tp, uptake by all Tp, T+, decay
            d4 = p[1] * Tpro[i] - mu[1,0] * Tpos[i] - mu[1,1] * Tpro[i] - lam[1] * test[i] + D[1] * (test[i - 1] + test[i + 1] - 2 * test[i])
            dtest.append(d4)
    
        #Equation for T+
        d0 = r[0] * Tpos[i] * (1 - ((Tpos[i] + Tpro[i] + Tneg[i])/(K + rho[0] * f_res(o2[i] , lim[0,0])*f_res(test[i],lim[1,0])))) - delta[0] * Tpos[i] + ovr(Tpos, i, K_m[0], M) 
           
        #Equation for Tp 
        d1 = r[1] * Tpro[i] * (1 - ((Tpos[i] + Tpro[i] + Tneg[i])/(K + rho[1] * f_res(o2[i], lim[0,1]) * f_res(test[i],lim[1,1])))) - delta[1] * Tpro[i] + ovr(Tpro, i, K_m[0], M)
                       
        #Equation for T-    
        d2 = r[2] * Tneg[i] * (1 - ((Tpos[i] + Tpro[i] + Tneg[i])/(K + rho[2] * f_res(o2[i],lim[0,2])))) - delta[2] * Tneg[i] + ovr(Tneg, i, K_m[0], M)
   
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
  
def ovr(cells, i, k_m, M):

    ovrflw = - max(cells[i] - k_m, 0)
    if i == 0:
        ovrflw += max( (cells[i + 1] - k_m)/2, 0)
    if i == M-1:
        ovrflw += max( (cells[i - 1] - k_m)/2, 0)
    if 0 < i < M-1:
        if i == 1:
            ovrflw += max( cells[i - 1] - k_m, 0) + max( (cells[i + 1] - k_m)/2, 0)
        if i == M-2:
            ovrflw += max( (cells[i - 1] - k_m)/2, 0) + max(cells[i + 1] - k_m, 0)
        else:    
            ovrflw += max( (cells[i - 1] - k_m)/2, 0) + max( (cells[i + 1] - k_m)/2, 0) 
    return ovrflw
 
