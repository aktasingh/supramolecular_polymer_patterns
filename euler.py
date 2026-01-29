#!/usr/bin/env python
# coding: utf-8

import numpy as np

print("Assigning parameters")

# Simulation parameters
tmax = 100  # total time
nstep =1000000  # total steps 
dt = tmax/nstep 
dx = 0.5  # lattice spacing 
x = np.arange(0, 100, dx)
t = np.arange(0, tmax, dt)
space = len(x) 
tim = len(t) 


# Model parameters
ka = 4  # monomer activation rate constant
kd2 = 5000  # polymer end deactivation rate constant
kd3 = 0.75  # polymer chain deactivation rate constant
kp = 1  # rate constant for assembly
ctot = 1e4  # total monomer concentration 


# Matrix to store variable data
dd = np.zeros((space, space))
aa1 = np.zeros((space, space))
mm1 = np.zeros((space, space))
mm0 = np.zeros((space, space))


# Initial condition 
m1 = np.zeros((space, space)) 
m0 = np.zeros((space, space))
diff_factor = np.zeros((space,space)) 

print("Setting up initial condition")

# (Random perturbation)
a1 = 5000 + 1e-2 * np.random.rand(space, space)
m1 = 86.0367 + 1e-2 * np.random.rand(space, space)
m0 = 1.96559 + 1e-2 * np.random.rand(space, space)
d = (ctot-a1-m1) + 1e-2 * np.random.rand(space, space) 

# (Gaussian perturbation)
a1 = 5000 + np.zeros((space, space))
xc = (space*dx)/2 
yc = (space*dx)/2 
inv_sigma = 0.5
for i in range(space): 
    for j in range (space):
        x = i*dx
        y = j*dx
        m1[i,j] = 0.086*np.exp(- ((x-xc)**2 + (y-yc)**2)/2*inv_sigma**2) + 86.0367
        m0[i,j] =0.0019*np.exp(- ((x-xc)**2 + (y-yc)**2)/2*inv_sigma**2) + 1.96559
d = (ctot-a1-m1) + np.zeros((space, space))


d1 = 1  # monomer diffusion coefficient 
diff_factor = m0/m1  # polymer diffusion coefficient

np.savez(f"step_0.0.npz",d=d,a1=a1,m1=m1,m0=m0) 


# Numerical integration using Euler method
print("Starting integration")
for k in range (1, tim):
    # deactivated monomer concentration
    dd[1:-1,1:-1] = ( d[1:-1,1:-1] + 
        dt * ((-ka * d[1:-1,1:-1]) + (2 * kd2 * m0[1:-1,1:-1])) +
        (dt/dx**2) * (d[0:-2,1:-1] + d[2:,1:-1] + d[1:-1,0:-2] + d[1:-1,2:] -4*d[1:-1, 1:-1])
        )
    # activated monomer concentration
    aa1[1:-1,1:-1] = (a1[1:-1,1:-1] +
        dt * ((ka * d[1:-1,1:-1]) - (2 * kp *a1[1:-1,1:-1] * m0[1:-1,1:-1])) +
        (dt/dx**2) * (a1[0:-2,1:-1] + a1[2:,1:-1] + a1[1:-1,0:-2] + a1[1:-1,2:] -4*a1[1:-1, 1:-1])
        )
    
    # polymer mass concentration
    mm1[1:-1,1:-1] = (m1[1:-1,1:-1] + 
        dt *((-2 * kd2 * m0[1:-1,1:-1]) + (2 * kp * a1[1:-1,1:-1] * m0[1:-1,1:-1])) +
        (dt/dx**2) * (((diff_factor[1:-1,1:-1]) * (m1[0:-2,1:-1] + m1[2:,1:-1] +m1[1:-1,0:-2] + m1[1:-1,2:] -4*m1[1:-1, 1:-1])) + 
        (0.25 * (((diff_factor[2:,1:-1]-diff_factor[0:-2,1:-1]) * (m1[2:,1:-1]-m1[0:-2,1:-1])) + ((diff_factor[1:-1,2:]-diff_factor[1:-1,0:-2]) * (m1[1:-1,2:]-m1[1:-1,0:-2])))))
        )
    
    # polymer number concentration
    mm0[1:-1,1:-1] = (m0[1:-1,1:-1] + 
        dt *((kd3 * (m1[1:-1,1:-1] - 4 * m0[1:-1,1:-1])) - (((2 * kd2 * 4.6 * m0[1:-1,1:-1]) / ((m1[1:-1,1:-1] / m0[1:-1,1:-1]) - 2)**2) + (2 * kd3 * m0[1:-1,1:-1]) + (kp * m0[1:-1,1:-1]**2))) +
        (dt/dx**2) * (((diff_factor[1:-1,1:-1]) * (m0[0:-2,1:-1] + m0[2:,1:-1] + m0[1:-1,0:-2] + m0[1:-1,2:] -4*m0[1:-1, 1:-1])) +
        (0.25 * (((diff_factor[2:,1:-1]-diff_factor[0:-2,1:-1]) * (m0[2:,1:-1]-m0[0:-2,1:-1])) + ((diff_factor[1:-1,2:]-diff_factor[1:-1,0:-2]) * (m0[1:-1,2:]-m0[1:-1,0:-2])))))
        )
           
    # No-flux boundary condition
    dd[0, 1:-1] = dd[1, 1:-1]
    dd[-1, 1:-1] = dd[-2, 1:-1]
    dd[1:-1, 0] = dd[1:-1, 1]
    dd[1:-1, -1] = dd[1:-1, -2]
    
    aa1[0, 1:-1] = aa1[1, 1:-1]
    aa1[-1, 1:-1] = aa1[-2, 1:-1]
    aa1[1:-1, 0] = aa1[1:-1, 1]
    aa1[1:-1, -1] = aa1[1:-1, -2]
    
    mm1[0, 1:-1] = mm1[1, 1:-1]
    mm1[-1, 1:-1] = mm1[-2, 1:-1]
    mm1[1:-1, 0] = mm1[1:-1, 1]
    mm1[1:-1, -1] = mm1[1:-1, -2]
    
    mm0[0, 1:-1] = mm0[1, 1:-1]
    mm0[-1, 1:-1] = mm0[-2, 1:-1]
    mm0[1:-1, 0] = mm0[1:-1, 1]
    mm0[1:-1, -1] = mm0[1:-1, -2]


    # Update boundary corners
    dd[0, 0], dd[0, -1], dd[-1, 0], dd[-1, -1] = dd[1, 1], dd[1, -2], dd[-2, 1], dd[-2, -2]
    aa1[0, 0], aa1[0, -1], aa1[-1, 0], aa1[-1, -1] = aa1[1, 1], aa1[1, -2], aa1[-2, 1], aa1[-2, -2]
    mm1[0, 0], mm1[0, -1], mm1[-1, 0], mm1[-1, -1] = mm1[1, 1], mm1[1, -2], mm1[-2, 1], mm1[-2, -2]
    mm0[0, 0], mm0[0, -1], mm0[-1, 0], mm0[-1, -1] = mm0[1, 1], mm0[1, -2], mm0[-2, 1], mm0[-2, -2]
    
    if np.remainder(k,1000)==0:
        print(f"Completed {k} steps", flush = True)
        np.savez(f"step_{k/10000}.npz",d=dd,a1=aa1,m1=mm1,m0=mm0) 
    
    d, a1, m1, m0 = dd.copy(), aa1.copy(), mm1.copy(), mm0.copy()
    
    diff_factor = m0/m1
    
print("Completed integration")



    

