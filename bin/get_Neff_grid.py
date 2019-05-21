import numpy as np
import matplotlib.pyplot as plt

for N_det in [1e4,1e5,1e6]:
    for fsky in [0.25,0.5,0.75]:
        for beam in [1,2,3,4]:
            fisher_matrix = np.loadtxt('../data/fisher_matrices/fisher_grid_{:d}_{:3.2f}_{:d}_p_l_nor.fish'.format(int(N_det),fsky,beam))[:-1,:-1]
            fisher_matrix[2,2] += 1/(0.01*69)**2
            #fisher_matrix[5,5] += 1e4
            #fisher_fixed = np.delete(fisher_matrix,2, axis = 0)
            #fisher = np.delete(fisher_fixed,2, axis = 1)
            cov = np.linalg.inv(fisher_matrix)
            print(np.sqrt(cov[-1,-1])*100, end = ' ')
        print('')