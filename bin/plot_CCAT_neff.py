import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

sensi_array = np.linspace(0.1,5,20)

error_list = []

def get_fisher_PyRayTE(fisher_file):
    fisher = np.loadtxt(fisher_file)
    fisher[5,5] += 1e4
    fisher_fixed = np.delete(fisher,(7),axis = (0))
    fisher_refixed = np.delete(fisher_fixed,(7),axis = (1))
    return fisher_refixed

fisher_nor = get_fisher_PyRayTE('../data/fisher_matrices/fisher_CCAT_v2_p_nol_nor.fish')
fisher_r = get_fisher_PyRayTE('../data/fisher_matrices/fisher_CCAT_v2_p_nol_r.fish')


error_nor = np.sqrt(np.linalg.inv(fisher_nor)[6,6])
error_r = np.sqrt(np.linalg.inv(fisher_r)[6,6])

print(error_nor,error_r)
print((error_nor-error_r)/error_nor *100)

for i in range(len(sensi_array)):
    file = '../data/fisher_matrices/fisher_CCAT_{:3.2f}_p_nol_nor.fish'.format(sensi_array[i])
    fisher = get_fisher_PyRayTE(file)
    error_list.append(np.sqrt(np.linalg.inv(fisher)[6,6]))

f = interp1d(np.array(error_list),sensi_array)
x_r = f(error_r)
x_nor = f(error_nor)
    
beta = x_nor / x_r
print(beta)
plt.rc('text',usetex = True)
plt.figure()
plt.plot(x_nor/sensi_array,error_list, c = 'k', lw = 2.5, label = r'$\sigma(N_\mathrm{eff})$ as a function of effort')
plt.plot([-1,x_nor/x_r],[error_r,error_r],c='r',lw = 2,ls = 'dashed', label = 'Forecast for CCAT + SO + PLANCK : TT+TE+EE + Rayleigh')
plt.plot([-1,x_nor/x_nor],[error_nor,error_nor],c='g',lw = 2,ls = 'dashed', label = 'Forecast for CCAT + SO + PLANCK : TT+TE+EE')
plt.plot([x_nor/x_r,x_nor/x_r],[0,error_r],c='r',lw = 2,ls = 'dashed')
plt.plot([x_nor/x_nor,x_nor/x_nor],[0,error_nor],c='g',lw = 2,ls = 'dashed')
plt.annotate(s = '', xy = (1,5e-2), xytext = (beta,5e-2), arrowprops = dict(arrowstyle = '<->'))
plt.annotate(s = '', xy = (0.55,error_r), xytext = (0.55,error_nor), arrowprops = dict(arrowstyle = '<->'))
plt.text(1.75,5.05e-2, r'$\times 2.6$',fontsize = 20)
plt.text(0.58,5.6e-2, r'$8\%$',fontsize = 20)
plt.yticks(list(plt.yticks()[0]) + [error_nor,error_r])
plt.ylim(4.5e-2,7e-2)
ax = plt.gca()
ax.ticklabel_format(axis = 'y', style = 'plain', useOffset = True, useMathText = True)
#plt.xscale('log')
plt.xlim(0.4,3.5)
plt.legend(prop = {'size' : 15})
plt.tick_params('both', labelsize = 18 )
plt.xlabel('Equivalent effort (eg. number of detectors, integration time, ...)', fontsize = 15)
plt.ylabel(r'$\sigma(N_\mathrm{eff})$', fontsize = 18)

plt.show()