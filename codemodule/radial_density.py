import numpy as np
import matplotlib.pyplot as plt
from .lennard_jones_v2 import d_matrix
import numba as nb
import codemodule as cm
from scipy import integrate
import itertools
def rad_frame(frame,dr,rho):
    #frame = frame.T
    print(frame[0])
    N_total = np.shape(frame)[1]
    L = (N_total/rho)**(1/3)
    d_arr = np.sqrt(d_matrix(frame,L))
    bins = int(L/(2*dr))
    #arb_indices = np.random.choice(N_total,size = n_part,replace = False)
    r_arr = np.array(list(d_arr[np.triu_indices(N_total,k=1)]))
    #counter = 0
    #for i in arb_indices:
     ##   for j in range(N_total-1):
         #   if i != j:
       #         r_arr[counter] = d_arr[i,j]
        #        counter += 1 
        
    hist, edges= np.histogram(r_arr,bins,range = (0,L/2))
    hist = 2*hist

    r = (edges[:-1]+edges[1:])/2
    norm = 4*np.pi*r**2*dr*rho*N_total
    g = hist/norm  

    
    return r,g

def radial_dist(xyz,dr,rho):
    #xyz as a list of all xyz frames
    g_list = []
    for i in range(len(xyz)):
        r,g = rad_frame(xyz[i],dr,rho)
        g_list.append(g)

    avg_g = np.mean(g_list,axis = 0)

    return r,avg_g

def radial_dist_types(frames,types,dr,rho):
    frames = np.transpose(frames)
    N_total = np.shape(frames)[1]
    L = (N_total/rho)**(1/3)
    d_arr = np.sqrt(d_matrix(frames,L))
    unique, counts = np.unique(types, return_counts = True)
    combinations = [(a, b) for idx, a in enumerate(unique) for b in unique[idx + 1:]]
    d_types = []
    for i in range(len(combinations)):
        arr = np.zeros((N_total,N_total))
        for j in range(N_total):
            for k in range(N_total):
                part_types = [types[i],types[j]]
                if part_types == combinations[i]:
                    arr[i,j] = d_arr[i,j]
        d_types.append(arr)

    print(d_types)

                    


    return combinations 
def pressure(g,r,rho,eps,sig,types,trunc):
    u_p = np.zeros(len(r))
    for i in range(len(r)):
        print(cm.lj_force(r[i],eps[0,0],sig[0,0])) #need to figure out a way to track types in r
        u_p[i] = cm.lj_force(r[i],eps[0,0],sig[0,0]) #need to figure out a way to track types in r

    integrand = r*u_p*g*4*np.pi*r**2
    r_space = np.linspace(0,2.5,num = len(r))
    P = rho-rho**2*integrate.simps(integrand,r_space) #idk how scipy works lmao or if i have it... gotta check when i have internet
    return P 






