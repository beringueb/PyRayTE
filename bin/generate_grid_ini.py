import numpy as np 

ini_file = '../ini/grid.ini'

for N_det in [1e4,1e5,1e6]:
    for fsky in [0.25,0.5,0.75]:
        for beam in [1,2,3,4]:
            with open(ini_file,'r') as f:
                text = f.read()
            new_text = text.replace('noise_file = ../data/noise_spectra/experiment_grid/noise_10000_0.75_1_lensing.dat','noise_file = ../data/noise_spectra/experiment_grid/noise_{:d}_{:3.2f}_{:d}_lensing.dat'.format(int(N_det),fsky,beam))
            new_new_text = new_text.replace('fsky = 0.5','fsky = {:3.2f}'.format(fsky))
            new_file = '../ini/grid/grid_{:d}_{:3.2f}_{:d}.ini'.format(int(N_det),fsky,beam)
            with open(new_file,'w') as f:
                f.write(new_new_text)
    
        
