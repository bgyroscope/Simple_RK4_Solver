# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 12:18:09 2017

Solving a system of coupled differential equations

@author: bgyroscope
"""

import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib.font_manager import FontProperties
    
#generate y_n+1
def nexty(fun, tn, yn, h):
    '''
    Produce the next point using RG4 integration, ie yn+1
    tn: current timestep
    yn: current array of yn values
    h: equals the time step size
    '''
    #define my points note that k's are arrays
    k1 = fun(tn, yn)
    k2 = fun(tn + h/2.0, yn + k1* h/2.0)
    k3 = fun(tn + h/2.0, yn + k2* h/2.0)
    k4 = fun(tn + h,     yn + k3 * h )
    
    return yn + h/6.0 * (k1 + 2*k2 + 2*k3 + k4)
    
def solve_coupled(f, n, interval, N, y0):
    """Solves a set of coupled ODEs Runge-Kutta method.

    f: set of equations f(t, y) where y is an array of variables
    n: integer number of equations 
    interval: array specifying where to solve, np.array([a,b])
    N: integer number of steps in the approximation
    y0: initial condition on x at t=a, given as an array
    

    Returns array of solution points
    """
    #defining time
    ti = interval[0]
    tf = interval[1]
    
    #step size
    h = 1.0*(tf- ti) / N
    
    time = np.linspace(ti, tf, N)
    
    sol = np.zeros([n, len(time) ]) #defining solution array
    
    sol[:,0] = y0  #setting initial conditions into the zeroth column
    
    #actual solving each step solves for n+1
    for j in range(N-1):
        yn = sol[:,j] #start with known yn
        tn = time[j]
        
        sol[:,j+1] = nexty(f, tn, yn, h)
        
    return sol


###########################Acutally starting this.~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define differential eqn.
def funExp(tn, yn):
    '''defines the dy/dt = f(t,y) '''
    return (yn)
    
    
#initial value problem. fun see above
y0 = 1;

#numerics
numpts = 200
start = 0
end = 5

t = np.linspace(start, end, numpts)
exact = np.exp(t)

#produce the numerical solution
interval = np.array([0,5])
y0 = np.array([1])
approx = solve_coupled( funExp, 1, interval, numpts, y0 )


plt.figure(1) 
plt.xlabel('t')
plt.grid(True)
# plt.hold(True)

lw = 1 #linewidth

plt.plot(t, exact, 'b', linewidth = lw)
plt.plot(t, approx[0], 'r', linewidth = lw)


plt.legend(('Exact', 'Approx'), prop=FontProperties(size=16))
plt.title('RG4 Method')

plt.show() 

