# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 17:46:37 2024

@author: emmac
"""
import numpy as np
from scipy import optimize

Ut_i_guess=1  #m/s, inital guess
tolerance=0.001  #convergence tolerance (%)
max_iter=50  #max iterations
g= 9.81 #m/s^2

#Design Column (client specs)
dp=1/1000  #m
rho_p=1650  #kg/m^3
rho_fluid=998.2  #kg/m^3,  online value for water at 20C
mew_fluid=0.0010016  #Pa*s, online value for water at 20C
V_flow=50/3600  #m^3/s

Diameter_guess=1  #m

def CD(Re):
    if Re<0.1:
        CD=24/Re
        return CD    #return calculated CD
    
    elif 0.1<Re<1000:
        CD=(24/Re)*(1+0.14*Re**0.70)
        return CD    #return calculated CD
    
    elif 1000<Re<3.5*10**5:
        CD=0.445
        return CD    #return calculated CD
    
    else:
        CD=0.445
        return CD    #return calculated CD

def Ut_calculator(Ut_guess): #calculate terminal velocity
    
    for i in range(max_iter):
        Re_guess=rho_fluid*Ut_guess*dp/mew_fluid  #using guessed Ut
        CD_guess=CD(Re_guess) #calculate drag coefficient
        Ut_new=((4*g*dp*(rho_p-rho_fluid))/(3*rho_fluid*CD_guess))**0.5
            
        if abs((Ut_guess-Ut_new)/Ut_new)*100<tolerance:  #check convergence (error in %)
            break #break out of for loop
        
        if i==max_iter:
            print("Maximum iterations reached, convergence not achieved")
            break #break out of for loop
        
        else:
            Ut_guess=Ut_new  #update Ut_guess if no convergence

    return CD_guess, Re_guess, Ut_guess   #converged CD, Re and Ut


def porosity(Dp, Dc):
    if Dp/Dc <= 0.265:
        E=0.4+(0.1*(np.exp(10.686*Dp/Dc)-1))
        return E    #return calculated porosity
    
    if 0.265<Dp/Dc<=0.538:
        E=0.846-(1.898*Dp/Dc+2.725*(Dp/Dc)**2)
        return E    #return calculated porosity
    
    if 0.538<Dp/Dc:
        E=1-((2/3)*((Dp/Dc)**3)/(np.sqrt((2*Dp/Dc)-1)))
        return E    #return calculated porosity

def solver_func(Dc):
    extension=0.285
    Eo=porosity(dp,Dc)
    E=1-((1-Eo)/(extension+1)) #porosity for column diameter of Dc input
    
    U_guess=V_flow/(np.pi*(Dc/2)**2)
    Ut_guess=Ut_calculator(Ut_i_guess)[2] #Converged terminal velocity w initial guess of 1m/s
    Re_guess=Ut_calculator(Ut_i_guess)[1]  #Converged Reynolds number w initial guess of 1m/s
    
    #n calculation based on Reynolds number
    if  Re_guess<0.2:
        n=4.65    #return calculated Richardson-Zaki index
    if 0.2<=Re_guess<1:
        n=4.4*Re_guess**(-0.03)    #return calculated Richardson-Zaki index
    if 1<=Re_guess<500:
        n=4.4*Re_guess**(-0.1)    #return calculated Richardson-Zaki index
    if Re_guess>=500:
        n=2.4    #return calculated Richardson-Zaki index
    
    #left and right side of richarson-zaki equation
    LS= E**n
    RS= U_guess/Ut_guess
    
    return LS-RS
    
#fsolve for LS and RS of Richardson-Zaki
Dc_guess=1
Dc_solved=optimize.fsolve(solver_func, Dc_guess, args=(), xtol=1e-8)

#Results
Eo=porosity(dp,Dc_solved[0])
E_final=1-((1-Eo)/(0.285+1))
U_final=V_flow/(np.pi*(Dc_solved[0]/2)**2)
Re_final=rho_fluid*U_final*dp/mew_fluid
n_final=4.4*Re_final**(-0.1)
CD_final=CD(Re_final)
Ut_final=(E_final**n_final)-U_final

print("The optimized diameter column is", Dc_solved[0], "m")
print("The final porosity is", E_final)
print("The final velocity is", U_final,"m/s")
print("The final Reynolds number is", Re_final)
print("The final Richardson-Zaki index is", n_final)
print("The final drag coefficient is", CD_final)
print("The final terminal velocity is", Ut_final,"m/s")
