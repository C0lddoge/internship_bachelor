import numpy as np



def write_xyz(f_out,xyz,types,counter,comment=''):  
    f_out.write('{:}\nsteps = {:},{:}\n'.format(np.shape(xyz)[1],counter,comment))
    for j in range(np.shape(xyz)[1]):
            f_out.write('{:} {:12.6f} {:12.6f} {:12.6}\n'.format(types[j],xyz[0,j],xyz[1,j],xyz[2,j]))
    pass


def write_txt(txt_file,energy,total_step):
    txt_file.write('{:}; {:12.6f}\n'.format(total_step,energy))
    pass

def write_xy(filepath,x,y):
    #takes x and y input arrays and creates a txt file 
    with open(filepath,'x') as f:
        for i in range(len(x)):
            x_ = x[i]
            y_ = y[i]
            f.write('{:}; {:}\n'.format(x_,y_))
    pass


def read_xyz(xyz_file,start_frame,end_frame):

    with open(xyz_file) as f_read:

        #N_total = int(f_read.readline().rstrip())
        #sframe_line = start_frame*N_total+2*start_frame
        #eframe_line = end_frame*N_total+2*end_frame
        frame_line = 69
        data = []
        coordinates = []

        for i, line in enumerate(f_read):
            if i%(frame_line) == 0:
                if i == 0:
                    frame_line = int(line.rstrip()) + 2
                    sframe_line = frame_line*start_frame
                    eframe_line = frame_line*(end_frame+1)
                elif i > sframe_line:
                    npcoordinates = np.array(coordinates)
                    data.append(npcoordinates)
                    coordinates.clear()
            elif i%frame_line != 1:
                if i > sframe_line:
                    t,x,y,z = line.split()
                    coordinates.append([float(x),float(y),float(z)])
                if i > eframe_line:

                    break
        types = np.loadtxt(xyz_file,skiprows = 2, max_rows = frame_line-2,usecols = 0, dtype = int)
        
    return data,types
                
            
def tail_correction_lj(r_c,N,rho,sig,eps):

    return rho/2 * 4*eps*(3*sig**12-44*np.pi*r_c**8*sig**6)/(33*r_c**11)

    

                






















