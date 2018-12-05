import numpy as np 

data_so = '/Users/benjamin_brng/Documents/Cambridge/PhD/Rayleigh/Data/Noise_spectra/noise_CCAT_lensing.dat'


noise = np.loadtxt(data_so,usecols = (0,2,3,4,10,11,12))
np.savetxt('../data/noise_spectra/template.dat',noise)