import numpy as np
import matplotlib.pyplot as plt
from .lennard_jones_v2 import d_matrix
import numba as nb
from .statistics import bin_ana_input, bin_ana
from scipy import integrate
def rad_frame(frame,dr,rho):
    #frame = frame.T
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
    hist = hist
    
    V = 4/3 * np.pi*(edges[1:]**3-edges[:-1]**3)
    r = (edges[:-1]+edges[1:])/2
    norm = V*rho*N_total
    g = hist/norm  
    
    return r,g

def g_r_error(g_series):
    #returns error array for each bin

    Nframe= len(g_series)
    Nbin = len(g_series[0])
    series_matrix = np.zeros((Nframe,Nbin))
    for i in range(Nframe):
        for j in range(Nbin):
            series_matrix[i,j] = g_series[i][j]

    rand_bin = np.random.randint(0,Nbin-1)
    Nb = bin_ana_input(series_matrix[:,rand_bin])
    error = []
    for i in range(Nbin):
        error.append(np.sqrt(bin_ana(series_matrix[:,i],Nb)))

    return error




def radial_dist(xyz,dr,rho):
    #xyz as a list of all xyz frames
    g_list = []
    for i in range(len(xyz)):
        r,g = rad_frame(xyz[i],dr,rho)
        g_list.append(g)
    
    error = g_r_error(g_list)
    avg_g = np.mean(g_list,axis = 0)

    return r,avg_g,error

def combinations(array):
    combinations = []
    for i in range(len(array)):
        for j in range(len(array)):
            combinations.append([array[i],array[j]])

    return combinations 
        

def radial_frames_types(frames,pairs,counts,types,dr,rho):
    frames = frames.T
    N_total = np.shape(frames)[1]
    L = (N_total/rho)**(1/3)
    bins = int(L/(2*dr))
    d_arr = np.sqrt(d_matrix(frames,L))
    N_pair = len(pairs)
    pair_counts = []
    for i in range(N_pair):
        pair_counts.append(counts[pairs[i][0]]*counts[pairs[i][1]])


    d_types = []
    for i in range(N_pair):
        arr = []
        for j in range(N_total):
            for k in range(N_total):
                part_types = [types[j],types[k]]
                if part_types == pairs[i] and j != k :
                    arr.append(d_arr[j,k])
        d_types.append(arr)
    
    g_pair = [] 
    r_pair = []

    for i in range(N_pair):
        hist,edges = np.histogram(d_types[i],bins, range = (0,L/2))
        r = (edges[:-1]+edges[1:])/2
        V = 4/3 * np.pi*(edges[1:]**3-edges[:-1]**3)
        norm = V*pair_counts[i]/L**3
        g = hist/norm
        g_pair.append(g)
        r_pair.append(r)


    return r_pair,g_pair


def rad_dist_types(frames,types,dr,rho):
    g_pair_frame = []
    unique, counts = np.unique(types, return_counts = True)
    pairs = combinations(unique)
    for i in range(len(frames)):
        r_pair,g_pair = radial_frames_types(frames[i],pairs,counts,types,dr,rho)
        g_pair_frame.append(g_pair)
    g_pair_avg =[]
    N_pair = len(pairs)
    error = []
    for j in range(N_pair):
        helper = []
        for k in range(len(frames)):
            helper.append(g_pair_frame[k][j])
        error.append(g_r_error(helper))
        g_pair_avg.append(np.mean(helper,axis=0))
    
    
    return pairs,r_pair,g_pair_avg,error



def pressure(g,r,rho,eps,sig,types,trunc):
    u_p = np.zeros(len(r))
    for i in range(len(r)):
        print(cm.lj_force(r[i],eps[0,0],sig[0,0])) #need to figure out a way to track types in r
        u_p[i] = cm.lj_force(r[i],eps[0,0],sig[0,0]) #need to figure out a way to track types in r

    integrand = r*u_p*g*4*np.pi*r**2
    r_space = np.linspace(0,2.5,num = len(r))
    P = rho-rho**2*integrate.simps(integrand,r_space) #idk how scipy works lmao or if i have it... gotta check when i have internet
    return P 






