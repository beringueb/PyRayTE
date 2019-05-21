import numpy as np 
import matplotlib.pyplot as plt 

freqs = [27,39,93,145,225,280,400,500]

param_list = ['ombh2','omch2','theta_MC','109A_s','n_s','tau','N_eff','mass_neutrinos']
   
def get_fisher_PyRayTE(fisher_file):
    fisher = np.loadtxt(fisher_file)
    #fisher[5,5] += 1e4
    fisher_fixed = np.delete(fisher,(7),axis = (0))
    fisher_refixed = np.delete(fisher_fixed,(7),axis = (1))
    return fisher_refixed
    
def get_fisher_full(fisher_file):
    fisher = np.loadtxt(fisher_file)
    fisher[5,5] += 1e4
    return fisher
    
for fsky in [0.1,0.25,0.5]:
    for v in ['v1','v2']:
        for str in ['r']:
            fisher_10000 = get_fisher_full('../data/fisher_matrices/fisher_CMB_HD_lensed_{:s}_7000_{:3.2f}_p_l_{:s}.fish'.format(v,fsky,str))
            i = 6
            param = param_list[i]
            error_10000 = np.sqrt(np.linalg.inv(fisher_10000)[i,i])
            if v == 'v1' : 
                noise = 0.25
            else : 
                noise = 0.5
            if str == 'nor':
                str_r = 'no rayleigh'
            else : 
                str_r = 'rayleigh'
    
            #print('{:3.2f} & {:3.2f} & {:s} & ${:4.2e}$ \\'.format(fsky,noise,str_r,error_10000))
            print('${:4.2e}$ \\'.format(error_10000))

lmax_list = np.logspace(np.log10(31),4,30).astype('int')
error_r = []
error_nor = []
error_PICO_r = []
error_PICO_nor = []
for lmax in lmax_list:
    fisher_r = get_fisher_full('../data/fisher_matrices/CMB_HD_lmax/fisher_CMB_HD_lensed_pol_v1_{:d}_p_l_r.fish'.format(lmax))
    fisher_nor = get_fisher_full('../data/fisher_matrices/CMB_HD_lmax/fisher_CMB_HD_lensed_pol_v1_{:d}_p_l_nor.fish'.format(lmax))
    fisher_PICO_r = get_fisher_full('../data/fisher_matrices/PICO_lmax/fisher_PICO_lensed_pol_v1_{:d}_p_l_r.fish'.format(lmax))
    fisher_PICO_nor = get_fisher_full('../data/fisher_matrices/PICO_lmax/fisher_PICO_lensed_pol_v1_{:d}_p_l_nor.fish'.format(lmax))
    i = 6
    param = param_list[i]
    error_r.append(np.sqrt(np.linalg.inv(fisher_r)[i,i]))
    error_nor.append(np.sqrt(np.linalg.inv(fisher_nor)[i,i]))
    error_PICO_r.append(np.sqrt(np.linalg.inv(fisher_PICO_r)[i,i]))
    error_PICO_nor.append(np.sqrt(np.linalg.inv(fisher_PICO_nor)[i,i]))
   
print(error_r[-1])   
print(error_nor[-1]) 
print(error_PICO_r[-1]) 
print(error_PICO_nor[-1]) 

plt.rc('text', usetex = True)
plt.figure()
plt.plot([10,11000],[0.03,0.03], alpha = 0.8, color = 'k', label = 'CMB-S4 ($\ell_{max} = 5000 $)')
plt.plot([10,11000],[0.00918,0.00918], alpha = 0.3, color = 'k', linestyle = 'dashed')
plt.plot([10,11000],[0.0145,0.0145], alpha = 0.3, color = 'k', linestyle = 'dashed')
plt.plot([10,11000],[0.0208,0.0208], alpha = 0.3, color = 'k', linestyle = 'dashed')
plt.plot([10,11000],[0.0373,0.0373], alpha = 0.3, color = 'k', linestyle = 'dashed')
plt.plot(lmax_list, error_r, label = 'CMB in HD, Rayleigh Scattering')
plt.plot(lmax_list, error_nor, label = 'CMB in HD')
plt.plot(lmax_list, error_PICO_r, label = 'PICO, Rayleigh Scattering')
plt.plot(lmax_list, error_PICO_nor, label = 'PICO')
plt.xlabel('$\ell_{max}$', fontsize = 26)
plt.ylabel('$\sigma (N_\mathrm{eff})$', fontsize = 26)
plt.xlim(600,10000)
plt.xscale('log')
plt.yscale('log')
plt.tick_params(axis='both', which='major', labelsize=16)

plt.yticks([0.00918,0.0145,0.0208,0.03,0.0373,0.05,0.1,0.15,0.3], [0.00918,0.0145,0.0208,0.03,0.0373,0.05,0.1,0.15,0.3])
plt.xticks([600,750,1000,1500,2500,5000,7500,10000],[600,750,1000,1500,2500,5000,7500,10000])

plt.legend(prop = {'size':17})
plt.show()

plt.show()

            









