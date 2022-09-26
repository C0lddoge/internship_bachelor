import numpy as np
import matplotlib.pyplot as plt
import random
import math
import os 
import numba as nb
from .data import *


@nb.jit(nopython = True)
def lj(r,eps,sig,E_c = 0):
    x = (sig/r)**6
    return 4*eps*((x**2)-x)-E_c

@nb.jit(nopython = True)
def lj_force(r,eps,sig):
    return -24*eps*(sig**6/r**7 - 2*sig**12/r**13)

def intialize(N,L):
    N_total = np.sum(N)
    N_prime = math.ceil(N_total**(1/3))
    a = L/N_prime

    counter = 0
    xyz = np.zeros((3,N_total))
    for i in range(N_prime):
        for j in range(N_prime):
            for k in range(N_prime):
                if counter < N_total:
                    xyz[0,counter] = a*(i+0.5)
                    xyz[1,counter] = a*(j+0.5)
                    xyz[2,counter] = a*(k+0.5)
                    counter += 1



    types = np.zeros(N_total, dtype = int)
    counter = 0
    for i in range(len(N)):
        for j in range(N[i]):
            types[counter] = i
            counter += 1


    return xyz,types



@nb.jit(nopython = True)
def d_finder(vec1,vec2,L):
    """finds distance of 2 vectors in the minimal image condition. Returns distances squared. Only take root once truncation has been accounted for."""
    d_squared = 0.0
    for i in range(3):
        d_i = vec1[i]-vec2[i]
        d_i -= L*round(d_i/L)

        d_squared += d_i**2
    
    return d_squared
@nb.jit(nopython = True)
def d_vec(vec1,vec2,L):
    d_vec = np.zeros(3)
    for i in range(3):
        d_vec[i] = vec1[i]-vec2[i]
        d_vec[i] -= L*round(d_vec[i]/L)

    return d_vec

@nb.jit(nopython = True)
def d_matrix(xyz,L):
    """ matrix that contains all distance squares for each particle"""
    N_total = np.shape(xyz)[1]
    d_matrix = np.zeros((N_total,N_total))
    for i in range(N_total-1):
        for j in range(i+1,N_total):
            r = d_finder(xyz[:,i],xyz[:,j],L)
            d_matrix[i,j] = r
            d_matrix[j,i] = r


    return d_matrix

def energy(types,p,xyz,trunc,sig,eps,N,L):
    N_total = np.shape(xyz)[1]
    E_shift = np.zeros(len(N))


    #call d_finder; dont use d_matrix
    for i in range(len(N)):
        if sig[types[p],i] != 0:
            E_shift[i] = lj(trunc*sig[types[p],i],eps[types[p],i],sig[types[p],i])
    
    energy_ = 0
    for i in range(N_total):
        if p != i:
            d_squared = d_finder(xyz[:,p],xyz[:,i],L)
            if d_squared<(trunc*sig[types[p],types[i]])**2:
                d = np.sqrt(d_squared)
                energy_ += lj(d,eps[types[p],types[i]],sig[types[p],types[i]],E_shift[types[i]])
            else:
                energy_ += 0
    return energy_

def total_energy(xyz,types,eps,sig,trunc,L):
    N_total = np.shape(xyz)[1]
    d_arr = d_matrix(xyz,L)
    E_shift = np.zeros(np.shape(sig))
    for i in range(np.shape(sig)[0]):
        for j in range(np.shape(sig)[1]):
            if sig[i,j] != 0:
                E_shift[i,j] = lj(trunc*sig[i,j],eps[i,j],sig[i,j])

    energy_total = 0
    for i in range(N_total):
        energy_ = 0
        for j in range(i):
            if i != j:
                if d_arr[i,j]<(trunc*sig[types[i],types[j]])**2:
                    d_arr[i,j] = np.sqrt(d_arr[i,j])
                    energy_ += lj(d_arr[i,j],eps[types[i],types[j]],sig[types[i],types[j]],E_shift[types[i],types[j]])
                else:
                    energy_ += 0
        energy_total += energy_
    return energy_total


def pressure(xyz,types,eps,sig,L,beta,trunc,rho):
    V = L**3
    N_total = np.shape(xyz)[1]
    a = rho/beta
    b = 1/(V*3)
    c = 0
    for i in range(N_total):
        for j in range(i+1,N_total):
            r = d_vec(xyz[:,i],xyz[:,j],L)
            r_2 = np.dot(r,r)
            if r_2 < trunc**2:

                F = r/np.sqrt(r_2) * lj_force(np.sqrt(r_2),eps[types[i],types[j]],sig[types[i],types[j]])
                c += np.dot(r,F)

    return a+b*c




def mc_move(xyz,types,trunc,sig,eps,beta,N,L,max_displacement):
    trial_index = np.random.randint(len(types))
    E1 = energy(types,trial_index,xyz,trunc,sig,eps,N,L)
    trial_xyz = np.copy(xyz)
    trial_xyz[:,trial_index] = trial_xyz[:,trial_index] + (np.random.random(3)-0.5)*max_displacement 
    E2 = energy(types,trial_index,trial_xyz,trunc,sig,eps,N,L)

    if E2<E1:
        return 1,trial_xyz

    else:
        boltzmann = np.exp(-beta*(E2-E1))
        sample = np.random.random()
        if boltzmann > sample:
            return 1,trial_xyz
        else:
            return 0,xyz

def simulate(N,eps,sig,beta,rho,mc_moves,trunc,max_displacement,e_step,w_step,p_step,xyz_switch,e_switch,p_switch,filename = '',energy_file = '',pressure_file = ''):
    N_total = np.sum(N)
    L = (N_total/rho)**(1/3)
    xyz,types = intialize(N,L)
    energy_steps = np.zeros(math.ceil(mc_moves/e_step))
    accept_count = 0
    reject_count = 0

    if xyz_switch:
        f_out = open(filename,'x')
    if e_switch:
        e_out = open(energy_file,'x')
    if p_switch : 
        p_out = open(pressure_file,'x')
    for i in range(mc_moves):
            print(i)
            yn,xyz = mc_move(xyz,types,trunc,sig,eps,beta,N,L,max_displacement)
            if yn:
                accept_count +=1
            else:
                reject_count +=1
            if xyz_switch:
                
                if i%w_step == 0:
                    write_xyz(f_out,xyz,types,i,L)
            if e_switch:
            
                if i%e_step== 0:
                    E_tot = total_energy(xyz,types,eps,sig,trunc,L)
                    write_txt(e_out,E_tot,i)
            if p_switch:
                
                if i%p_step == 0:
                    p = pressure(xyz,types,eps,sig,L,beta,trunc,rho)
                    write_txt(p_out,p,i)

    if xyz_switch:

        f_out.close()

    if e_switch:

        e_out.close()
    if p_switch:

        p_out.close()
                
    return  accept_count/mc_moves



