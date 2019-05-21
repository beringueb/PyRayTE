import numpy as np 

ini_file = '../ini/CCAT_init.ini'

for sensi in np.linspace(0.1,5,20):
    with open(ini_file,'r') as f:
        text = f.read()
    new_text = text.replace('noise_file = ../data/noise_spectra/CCAT_v2.dat','noise_file = ../data/noise_spectra/CCAT_{:3.2f}.dat'.format(sensi))
    new_file = '../ini/CCAT/CCAT_{:3.2f}.ini'.format(sensi)
    with open(new_file,'w') as f:
        f.write(new_text)
        
        
