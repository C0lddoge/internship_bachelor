import numpy as np

# hard-coded constans - don't touch
L_box = 10 # box length
disk_d = 1.0 # disk diameter

def mc_move_dynamic(x, L_box, N_part, disk_d, max_displacement):
  '''
  MC move, where you add a small random displacement to a single particle to generate a trial configuation
  '''
  trial_x = np.copy(x) # copy of the coordinates
  trial_index = np.random.randint(N_part) # which disk to move
  trial_x[trial_index] += (np.random.random((2)) - 0.5) * max_displacement # add displacment

  # if there is an overlap, reject the trial, else return the trial
  if check_for_overlaps(trial_x, L_box, N_part, disk_d):
    return 0, x
  else:
    return 1, trial_x % L_box # wrap the particles through the PBC to the original box

def mc_move_static(x, L_box, N_part, disk_d, max_displacement = None):
  '''
  naive MC move, where you just place all of the particles at random
  '''
  trial_x = np.random.random((N_part, 2)) * L_box # initialize the trial configuration
  # if there is an overlap, reject the trial, else return the trial
  if check_for_overlaps(trial_x, L_box, N_part, disk_d):
    return 0, x
  else:
    return 1, trial_x

def check_for_overlaps(x, L_box, N_part, disk_d):
  '''
  return True if there is a pair of overlapping rings, or False if there are no overlaps
  '''
  # loop over all pairs of disks
  for p_id1 in range(N_part-1):
    for p_id2 in range(p_id1+1, N_part):
      d2 = 0.0 # buffer for distance squared
      # loop over dimensions
      for d in range(2):
        dxk = x[p_id1, d] - x[p_id2, d] # displacement in the dimension d
        dxk -= L_box * round(dxk / L_box) # apply the periodic boundary conditions
        d2 += dxk * dxk # sum over dimensions
      # if there is a single overlap, where are done
      if d2 < disk_d**2:
        return True
  return False 

def write_xyz(f_out, x, N_part, comment = ''):
  '''
  write one xyz frame to an open file.
  '''
  f_out.write('{:d}\n{:}\n'.format(N_part, comment)) # comment line
  # one line per atom, z coordinate is fixed as 0
  for i in range(N_part):
    f_out.write('X {:12.6f} {:12.6f} 0.0\n'.format(*x[i]))