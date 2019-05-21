import numpy as np
import matplotlib.pyplot as plt 

param_list = ['ombh2','omch2','theta_MC','109A_s','n_s','tau','N_eff','mass_neutrinos']
param_list_joel = ['tau','omch2','109A_s','theta_MC','ombh2','n_s']

file_joel = "../data/templates/LCDM_V3.txt"


def read_error_joel(file):
    error_dict = {}
    for i in range(9):
        error = np.loadtxt(file)[i,1:]
        for param in param_list_joel:
            if i == 0:
                error_dict[param] = []
            index = param_list_joel.index(param)
            error_dict[param].append(error[index])
    return error_dict
    
    
def get_fisher_PyRayTE(fisher_file):
    fisher = np.loadtxt(fisher_file)
    fisher_fixed = np.delete(fisher,(6,7),axis = (0))
    fisher_refixed = np.delete(fisher_fixed,(6,7),axis = (1))
    return fisher_refixed
    
def get_error_ben(error_dict,fisher,i):
    cov = np.linalg.inv(fisher)

    for param in param_list_joel:
        if i == 0:
            error_dict[param] = []
        index = param_list.index(param)
        error_dict[param].append(np.sqrt(cov[index,index]))
    return error_dict
    
i = 0
error_dict_ben = {}
for sensi in [0,1,2]:
    for fsky in [0.1,0.2,0.4]:
        fisher1 = get_fisher_PyRayTE('../data/fisher_matrices/fisher_SO+PLANCK_{:d}_{:3.1f}_p_l_nor.fish'.format(sensi,fsky))
        #fisher2= get_fisher_PyRayTE('../data/fisher_matrices/fisher_PLANCK_{:d}_{:3.1f}_p_l_nor.fish'.format(sensi,fsky))
        fisher = fisher1#+fisher2
        fisher[5,5] += 1e4
        error_dict_ben = get_error_ben(error_dict_ben, fisher, i)
        i += 1
        
error_dict_joel = read_error_joel(file_joel)
for value in range(len(error_dict_joel['109A_s'])):
    error_dict_joel['109A_s'][value] *= 1E9

fisher_planck_only = get_fisher_PyRayTE('../data/fisher_matrices/fisher_PLANCK_only_p_l_nor.fish')

cov_planck = np.linalg.inv(fisher_planck_only)
for param in param_list[:-2]:
    i = param_list.index(param)
    #print('Error on {} : {:3.1e}'.format(param,np.sqrt(cov_planck[i,i])))

for param in param_list_joel:
    print('Error on {:s} : {:4.2f} %'.format(param, ((np.asarray(error_dict_joel[param]) - np.asarray(error_dict_ben[param]))/np.asarray(error_dict_joel[param])*100).mean()))
plt.rc('text', usetex = True)
f,axs = plt.subplots(2,3, sharex = True)
line = 0
column = 0
for param in param_list_joel:
    axs[line,column].set_title(r'${}$'.format(param))
    axs[line,column].plot([0.1,0.2,0.4],error_dict_joel[param][0:3],c = 'r')
    axs[line,column].plot([0.1,0.2,0.4],error_dict_ben[param][0:3],c = 'r', ls = 'dashed')
    axs[line,column].plot([0.1,0.2,0.4],error_dict_joel[param][3:6],c = 'b')
    axs[line,column].plot([0.1,0.2,0.4],error_dict_ben[param][3:6],c = 'b', ls = 'dashed')
    axs[line,column].plot([0.1,0.2,0.4],error_dict_joel[param][6:],c = 'g')
    axs[line,column].plot([0.1,0.2,0.4],error_dict_ben[param][6:],c = 'g', ls = 'dashed')
    column += 1
    if column == 3 : 
        line = 1
        column = 0 
plt.show()
            
        