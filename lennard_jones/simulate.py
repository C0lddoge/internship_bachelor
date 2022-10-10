import numpy as np
import os
import codemodule as cm
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument(
            "--density",
            type = float,
            required = True
            )
parser.add_argument(
            "--N",
            nargs = '+',
            type = int 
            )

parser.add_argument(
        "--T",
        type = float,
        default = 1.0
        )

parser.add_argument(
        "--xyz_switch",
        type = bool,
        default = False 

        )
parser.add_argument(
        "--e_switch",
        type = bool,
        default = False

        )
parser.add_argument(
        "--p_switch",
        type = bool,
        default = False
        )

input_dictionary = vars(parser.parse_args())

def sim_name(N, density,T):
    print(N)
    N_total = sum(N)
    types = len(N)
    return 'LJ_N{:d}_types_{:d}_dens{:.2e}_T{:.2e}'.format(N_total,types,density, T)

def simulate(N, density,T,xyz_switch,e_switch,p_switch):
    if not os.path.isdir('data'):
      os.mkdir('data')
    sim_path = '{:}/data/{:}/'.format(os.getcwd(), sim_name(N, density,T))
    if not os.path.exists(sim_path):
        os.mkdir(sim_path)
    plot_path = os.path.join(sim_path,'plots')
    if not os.path.isdir(plot_path):
        os.mkdir(plot_path)
    if not os.path.isdir(sim_path + 'results'):
        os.mkdir(sim_path+'results')
    E_tot_file = sim_path + '/ener.dat'
    xyz_file = sim_path + '/frames.xyz'
    pressure_file = sim_path +'/pressure.dat'

    beta = 1.0/T
    eps = np.array([[1,0],[0,1]])
    sig = np.array([[1,0],[0,1]])
    trunc = 2.5
    mc_moves = 10000
    max_displacement = 0.5
    energy_steps = 10
    w_steps = 10
    p_steps = 10
    with open(sim_path + '/log.dat', 'w') as out_f:
        out_f.write('\ntrunc = {:}\nmc_moves = {:}\nenergy_steps = {:}\n xyz_steps = {:}\npressure_steps = {:}\nN = {:}\ndensity = {:}\nT = {:}\n'.format(trunc,mc_moves,energy_steps,w_steps,p_steps,N, density,T))  
    out_f.close()
    cm.simulate(N,eps,sig,beta,density,mc_moves,trunc,max_displacement,energy_steps,w_steps,p_steps,xyz_switch,e_switch,p_switch,xyz_file,E_tot_file,pressure_file)


simulate(**input_dictionary)


















