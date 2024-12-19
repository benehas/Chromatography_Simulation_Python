# -*- coding: utf-8 -*-
"""
Created on Sat Aug 15 18:32:32 2020

@author: Benedict
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Aug 15 17:15:14 2020

@author: Benedict
"""

import numpy as np
import matplotlib.pyplot as plt
import copy as cp
import time 
#column Parameters
L = 15 #cm
tmax = 30 #s
u_inter = 1.5 #cm/s
epsilon = 0.72 
t_slug=2
keff=0.01
Rp=0.0025
Lambda=2

const_1 = (1-epsilon)/epsilon#constants to reduce computations
const_2 = 3/Rp #constants to reduce computations
D_ax = 0.05 #cm^2/s

#Numeric Parameters
dt = 0.005 #s
dx = 0.1 #cm
N=int(L/dx)
M=int(tmax/dt)


#Components
isotherme='Langmuir'
C=2 # Number of Components
params=np.ndarray(C,dict)

if isotherme=='Langmuir':
    params[0]={'KL': 1,'qmax': 0.5,'c_in':0.2}
    params[1]={'KL': 0.1,'qmax': 0.5,'c_in':0.2}

if isotherme=='Henry':
    params[0]={'Kh': 1,'c_in':1}
    params[1]={'Kh': 0.1,'c_in':1}

if isotherme=='SMA':
    params[0]={'Ksma': 3,'sigma':2,'ny':2,'c_in':2}
    params[1]={'Ksma': 1,'sigma':4,'ny':3,'c_in':2}

#Memory
concentration=np.zeros((C,N))
q=np.zeros((C,N))
if isotherme=='SMA':
    c_salt=np.zeros(N)
    c_salt_save = np.zeros(N)
concentration_s=np.zeros((C,N))
q_s=np.zeros((C,N))

#visualiation
plot=np.zeros((C,M))
t=[i*dt for i in range(M)]
def compute_dc_dx(c_i,i):
    return (c_i[i+1]-c_i[i-1])/(2*dx)

def compute_dc2_dx2(c_i,i):
    return (c_i[i+1]+c_i[i-1]-2*c_i[i])/(dx*dx)

def compute_cp_LG(c,q,i,params):
    A=np.ndarray((C,C))
    D=np.ndarray(C)
    
    for k in range(C):
        for l in range(C):
            if k==l:
                A[k][l] = params[l]['KL']*(q[k][i]-params[l]['qmax'])
            else:
                A[k][l]=q[k][i]*params[l]['KL']
        D[k]=-q[k][i]
    return D.dot(np.linalg.inv(A))

def compute_cp_SMA(q,k,i,params,csalt):
    tmp=0
    for r in range(C):
        tmp+=q[r][i]*(params[r]['ny']+params[r]['sigma'])
    return (q[k][i]/params[k]['Ksma'])*(csalt/(Lambda-tmp))**params[k]['ny']

#Main Computation

if __name__=='__main__':
    start=time.time()
    #Initial conditions
    for k in range(C):
        concentration[k][0]=params[k]['c_in']
    if isotherme=='SMA':
        c_salt[0]=1
        
    #Main Looper
    for i in range(M):
       concentration_s=cp.deepcopy(concentration) 
       q_s=cp.deepcopy(q)
       
       if isotherme=='SMA':
           c_salt_save = cp.deepcopy(c_salt)
       
       if isotherme=='Langmuir':
           for l in range(1,N-1):
               tmp=compute_cp_LG(concentration_s,q_s,l,params)
               for k in range(C):
                   concentration[k][l]+=(D_ax*compute_dc2_dx2(concentration_s[k],l)-u_inter*compute_dc_dx(concentration_s[k],l)-const_1*const_2*keff*(concentration_s[k][l]-tmp[k]))*dt
                   q[k][l]+=const_1*const_2*keff*(concentration_s[k][l]-tmp[k])*dt

       if isotherme=='Henry':
           for l in range(1,N-1):
               for k in range(C):
                   concentration[k][l]+=(D_ax*compute_dc2_dx2(concentration_s[k],l)-u_inter*compute_dc_dx(concentration_s[k],l)-const_1*const_2*keff*(concentration_s[k][l]-q[k][l]/params[k]['Kh']))*dt
                   q[k][l]+=const_1*const_2*keff*(concentration_s[k][l]-q[k][l]/params[k]['Kh'])*dt

       if isotherme=='SMA':
           for l in range(1,N-1):
               c_salt[l]+=(D_ax*compute_dc2_dx2(c_salt_save,l)-u_inter*compute_dc_dx(c_salt_save,l))*dt
               for k in range(C):
                   concentration[k][l]+=(D_ax*compute_dc2_dx2(concentration_s[k],l)-u_inter*compute_dc_dx(concentration_s[k],l)-const_1*const_2*keff*(concentration_s[k][l]-compute_cp_SMA(q,k,l,params,c_salt[l])))*dt
                   q[k][l]+=const_1*const_2*keff*(concentration_s[k][l]-compute_cp_SMA(q,k,l,params,c_salt[l]))*dt

        #boundary conditions
       for k in range(C):
           if i*dt<t_slug:
               concentration[k][0]=params[k]['c_in']
           else:
               concentration[k][0]=0
           concentration[k][-1]=concentration[k][-2]
       if isotherme=='SMA':
           c_salt[-1]=c_salt[-2]
           
        #add detecotor output for visualization
       for k in range(C):   
           plot[k][i]=concentration[k][-2]
       if 100*i/M%10==0:
        print('Simulation progress: '+str(100*i/M)+'%')
        print('Simulation time: '+str(i*dt)+'s')
        print('Phys time: '+str(time.time()-start)+'s')
    
    print('Total Computaion time: '+str(time.time()-start))
    #visualization
    plt.xlabel('time in s')
    plt.ylabel('C')
    for k in range(C):
        plt.plot(t,plot[k,:])
    plt.show()
        
          