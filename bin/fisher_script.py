import numpy as np 
import matplotlib.pyplot as plt 

freqs = [270,350,405]
beam = [0.8,0.6,0.5]
noise_v1 = [13,71,297]
noise_v2 = [7.3,41,172]

file_root = '../data/power_spectra/CCAT_v2+SO+PLANCK+PLANCK_low/fiducial_scalCls_lensed_'
ell,cmb = np.loadtxt(file_root + '1_1', usecols = (0,1), unpack = True)
ell,cross_270 = np.loadtxt(file_root + '1_5', usecols = (0,1), unpack = True)
ell,cross_350 = np.loadtxt(file_root + '1_7', usecols = (0,1), unpack = True)
ell,cross_405 = np.loadtxt(file_root + '1_8', usecols = (0,1), unpack = True)

ell,auto_270 = np.loadtxt(file_root + '5_5', usecols = (0,1), unpack = True)
ell,auto_350 = np.loadtxt(file_root + '7_7', usecols = (0,1), unpack = True)
ell,auto_405 = np.loadtxt(file_root + '8_8', usecols = (0,1), unpack = True)

noise_prim = ell*(ell+1)/2/np.pi * (5.*np.pi/(60*180))**2 * np.exp(ell*(ell+1) * (1. * np.pi / 10800.) ** 2 /(8. * np.log(2)))

noise_v1_270 = ell*(ell+1)/2/np.pi * (noise_v1[0]*np.pi/(60*180))**2 * np.exp(ell*(ell+1) * (beam[0] * np.pi / 10800.) ** 2 /(8. * np.log(2)))
noise_v1_350 = ell*(ell+1)/2/np.pi * (noise_v1[1]*np.pi/(60*180))**2 * np.exp(ell*(ell+1) * (beam[1] * np.pi / 10800.) ** 2 /(8. * np.log(2)))
noise_v1_405 = ell*(ell+1)/2/np.pi * (noise_v1[2]*np.pi/(60*180))**2 * np.exp(ell*(ell+1) * (beam[2] * np.pi / 10800.) ** 2 /(8. * np.log(2)))
noise_v2_270 = ell*(ell+1)/2/np.pi * (noise_v2[0]*np.pi/(60*180))**2 * np.exp(ell*(ell+1) * (beam[0] * np.pi / 10800.) ** 2 /(8. * np.log(2)))
noise_v2_350 = ell*(ell+1)/2/np.pi * (noise_v2[1]*np.pi/(60*180))**2 * np.exp(ell*(ell+1) * (beam[1] * np.pi / 10800.) ** 2 /(8. * np.log(2)))
noise_v2_405 = ell*(ell+1)/2/np.pi * (noise_v2[2]*np.pi/(60*180))**2 * np.exp(ell*(ell+1) * (beam[2] * np.pi / 10800.) ** 2 /(8. * np.log(2)))


error_270 = np.sqrt((cross_270**2 + (cmb + noise_prim)*(auto_270 + noise_v1_270) ) / 0.21/(2*ell+1.))
error_350 = np.sqrt((cross_350**2 + (cmb + noise_prim)*(auto_350 + noise_v1_350) ) / 0.21/(2*ell+1.))
error_405 = np.sqrt((cross_405**2 + (cmb + noise_prim)*(auto_405 + noise_v1_405) ) / 0.21/(2*ell+1.))

error_v2_270 = np.sqrt((cross_270**2 + (cmb + noise_prim)*(auto_270 + noise_v2_270) ) / 0.21/(2*ell+1.))
error_v2_350 = np.sqrt((cross_350**2 + (cmb + noise_prim)*(auto_350 + noise_v2_350) ) / 0.21/(2*ell+1.))
error_v2_405 = np.sqrt((cross_405**2 + (cmb + noise_prim)*(auto_405 + noise_v2_405) ) / 0.21/(2*ell+1.))

ell_bin = np.arange(0,5000,50).astype('int')
f,axs = plt.subplots(3,1, sharex = True)
plt.rc('text', usetex = True) 
axs[0].plot(ell,cross_270, color = 'r')
eb1 = axs[0].errorbar(x = ell_bin,y = cross_270[ell_bin], yerr = error_270[ell_bin]/np.sqrt(50), ls = 'None', marker = '+', color = 'r' ,capsize = 2, label = '270GHz, option 1')
eb2 = axs[0].errorbar(x = ell_bin,y = cross_270[ell_bin], yerr = error_v2_270[ell_bin]/np.sqrt(50), ls = 'None', marker = '+', color = 'r',capsize = 2, label = '270GHz, option 2')
eb1[-1][0].set_linestyle('dotted')
axs[0].set_xlim(0,3000)
axs[0].legend(loc = 'lower right', prop = {'size' : 14})

axs[1].plot(ell,cross_350, color = 'g')
eb1 = axs[1].errorbar(x = ell_bin,y = cross_350[ell_bin], yerr = error_350[ell_bin]/np.sqrt(50), ls = 'None', marker = '+', color = 'g'  ,capsize = 2, label = '350GHz, option 1')
eb2 = axs[1].errorbar(x = ell_bin,y = cross_350[ell_bin], yerr = error_v2_350[ell_bin]/np.sqrt(50), ls = 'None', marker = '+', color = 'g',capsize = 2, label = '350GHz, option 2')
eb1[-1][0].set_linestyle('dotted')
axs[1].set_xlim(0,3000)
axs[1].legend(loc = 'lower right', prop = {'size' : 14})

axs[2].plot(ell,cross_405, color = 'b')
ell_bin = np.arange(0,5000,50).astype('int')
eb1 = axs[2].errorbar(x = ell_bin,y = cross_405[ell_bin], yerr = error_405[ell_bin]/np.sqrt(50), ls = 'None', marker = '+', color = 'b' ,capsize = 2, label = '405GHz, option 1')
eb2 = axs[2].errorbar(x = ell_bin,y = cross_405[ell_bin], yerr = error_v2_405[ell_bin]/np.sqrt(50), ls = 'None', marker = '+', color = 'b',capsize = 2, label = '405GHz, option 2')
eb1[-1][0].set_linestyle('dotted')
axs[2].set_xlim(0,3000)
axs[2].legend(loc = 'lower right', prop = {'size' : 14})
axs[1].set_ylabel('$\ell*(\ell+1)/2\pi C_\ell^{Cross}$', fontsize = 16)
#axs[2].set_ylim(-100,5)
axs[2].set_xlabel('$\ell$', fontsize = 16)

plt.figure()
plt.plot(ell,cross_270, color = 'r')
eb1 = plt.errorbar(x = ell_bin,y = cross_270[ell_bin], yerr = error_270[ell_bin]/np.sqrt(50), ls = 'None', marker = '+', color = 'r' ,capsize = 2, label = '270GHz, option 1')
eb2 = plt.errorbar(x = ell_bin,y = cross_270[ell_bin], yerr = error_v2_270[ell_bin]/np.sqrt(50), ls = 'None', marker = '+', color = 'r',capsize = 2, label = '270GHz, option 2')
eb1[-1][0].set_linestyle('dotted')
plt.xlim(0,3000)
plt.legend(loc = 'lower right', prop = {'size' : 18})
plt.xlabel('$\ell$', fontsize = 25)
plt.ylabel('$\ell*(\ell+1)/2\pi C_\ell^{Cross}$', fontsize = 25)


plt.figure()
ell_bin = np.arange(0,5000,100).astype('int')

plt.plot(ell,cross_270, color = 'r',alpha = 0.75)
#eb1 = plt.errorbar(x = ell_bin,y = cross_270[ell_bin], yerr = error_270[ell_bin]/np.sqrt(200), ls = 'None', marker = '+', color = 'r' ,capsize = 2, label = '270GHz, option 1')
eb2 = plt.errorbar(x = ell_bin,y = cross_270[ell_bin], yerr = error_v2_270[ell_bin]/np.sqrt(200), ls = 'None', marker = '+', color = 'r',capsize = 2, label = '270GHz')
#eb1[-1][0].set_linestyle('dotted')
plt.plot(ell,cross_350, color = 'g',alpha = 0.75)
#eb1 = plt.errorbar(x = ell_bin,y = cross_350[ell_bin], yerr = error_350[ell_bin]/np.sqrt(200), ls = 'None', marker = '+', color = 'g' ,capsize = 2, label = '350GHz, option 1')
eb2 = plt.errorbar(x = ell_bin,y = cross_350[ell_bin], yerr = error_v2_350[ell_bin]/np.sqrt(200), ls = 'None', marker = '+', color = 'g',capsize = 2, label = '350GHz')
#eb1[-1][0].set_linestyle('dotted')
plt.xlim(0,3000)
plt.legend(loc = 'lower right', prop = {'size' : 18})
plt.xlabel('$\ell$', fontsize = 25)
plt.ylabel('$\ell*(\ell+1)/2\pi C_\ell^{Cross}$', fontsize = 25)


f,axs = plt.subplots(2,1) 
#axs[0].plot(ell, np.abs(cross_270/error_270), label = '270GHz v1')
#axs[1].plot(ell[500:], np.cumsum(np.abs(cross_270[500:]/error_270[500:])))
#axs[0].plot(ell, np.abs(cross_350/error_350), label = '350Ghz v1')
#axs[1].plot(ell[500:], np.cumsum(np.abs(cross_350[500:]/error_350[500:])))
axs[0].plot(ell, np.abs(cross_270/error_v2_270), label = '270GHz',c = 'r')
axs[1].plot(ell[500:], np.cumsum(np.abs(cross_270[500:]/error_v2_270[500:])),c='r')
axs[0].plot(ell, np.abs(cross_350/error_v2_350), label = '350Ghz',c='g')
axs[1].plot(ell[500:], np.cumsum(np.abs(cross_350[500:]/error_v2_350[500:])),c = 'g')
axs[1].set_xlabel('$\ell$', fontsize = 22)
axs[1].set_ylabel('Cumulative S/N', fontsize = 22)
axs[0].set_ylabel('S/N', fontsize = 22)
axs[0].legend(loc = 'upper right', prop = {'size':18})

axs[0].set_xlim(0,3000)
axs[1].set_xlim(500,3000)




plt.show()

param_list = ['ombh2','omch2','theta_MC','109A_s','n_s','tau','N_eff','mass_neutrinos']
   
def get_fisher_PyRayTE(fisher_file):
    fisher = np.loadtxt(fisher_file)
    fisher[5,5] += 1e4
    fisher_fixed = np.delete(fisher,(6),axis = (0))
    fisher_refixed = np.delete(fisher_fixed,(6),axis = (1))
    return fisher_refixed
    
def get_fisher_full(fisher_file):
    fisher = np.loadtxt(fisher_file)
    fisher[5,5] += 1e4
    return fisher

fisher_v1_opt1 = get_fisher_PyRayTE('../data/fisher_matrices/fisher_CCAT_v1_opt1_p_nol_r.fish')
fisher_v1_opt2 = get_fisher_PyRayTE('../data/fisher_matrices/fisher_CCAT_v1_opt2_p_nol_r.fish')
fisher_v2 = get_fisher_PyRayTE('../data/fisher_matrices/fisher_CCAT_v2_p_nol_r.fish')
fisher_SO = get_fisher_PyRayTE('../data/fisher_matrices/fisher_SO_p_nol_r.fish')
fisher_SO_nor = get_fisher_PyRayTE('../data/fisher_matrices/fisher_SO_p_nol_nor.fish')

for i in range(6):
    param = param_list[i]
    error_nor = np.sqrt(np.linalg.inv(fisher_SO_nor)[i,i])
    error_SO_r = np.sqrt(np.linalg.inv(fisher_SO)[i,i])
    error_v1_o1_r = np.sqrt(np.linalg.inv(fisher_v1_opt1)[i,i])
    error_v1_o2_r = np.sqrt(np.linalg.inv(fisher_v1_opt2)[i,i])
    error_v2_r = np.sqrt(np.linalg.inv(fisher_v2)[i,i])
    
    impvt1 = -(error_SO_r - error_nor)/error_nor*100
    impvt2 = -(error_v1_o1_r - error_SO_r)/error_SO_r*100
    impvt3 = -(error_v1_o2_r - error_SO_r)/error_SO_r*100
    impvt4 = -(error_v2_r - error_SO_r)/error_SO_r*100
    #print('Error on {:s}, SO nor :{:4.2e}, SO {:4.2e} ({:3.1f}%), v1 opt1 {:4.2e} ({:3.1f}%), v1 opt2 {:4.2e} ({:3.1f}%), v2 {:4.2e} ({:3.1f}%) \n'.format(param,error_nor,error_SO_r,impvt1,error_v1_o1_r,impvt2,error_v1_o2_r,impvt3,error_v2_r,impvt4))
    #print('Error on {:s}, SO nor :${:4.2e}$, SO ${:4.2e}$ (${:3.1f}%$), v1 opt1 ${:4.2e}$ (${:3.1f}%$), v1 opt2 ${:4.2e}$ (${:3.1f}%$), v2 ${:4.2e}$ (${:3.1f}%$) \n'.format(param,error_nor,error_SO_r,impvt1,error_v1_o1_r,impvt2,error_v1_o2_r,impvt3,error_v2_r,impvt4))
    print('${:4.2e}$ &'.format(error_v2_r), end = ' ')

fisher_v1_opt1 = get_fisher_full('../data/fisher_matrices/fisher_CCAT_v1_opt1_p_nol_r.fish')
fisher_v1_opt2 = get_fisher_full('../data/fisher_matrices/fisher_CCAT_v1_opt2_p_nol_r.fish')
fisher_v2 = get_fisher_full('../data/fisher_matrices/fisher_CCAT_v2_p_nol_r.fish')
fisher_SO = get_fisher_full('../data/fisher_matrices/fisher_SO_p_nol_r.fish')
fisher_SO_nor = get_fisher_full('../data/fisher_matrices/fisher_SO_p_nol_nor.fish')

i = 6
param = param_list[i]
error_nor = np.sqrt(np.linalg.inv(fisher_SO_nor)[i,i])
error_SO_r = np.sqrt(np.linalg.inv(fisher_SO)[i,i])
error_v1_o1_r = np.sqrt(np.linalg.inv(fisher_v1_opt1)[i,i])
error_v1_o2_r = np.sqrt(np.linalg.inv(fisher_v1_opt2)[i,i])
error_v2_r = np.sqrt(np.linalg.inv(fisher_v2)[i,i])

impvt1 = -(error_SO_r - error_nor)/error_nor*100
impvt2 = -(error_v1_o1_r - error_SO_r)/error_SO_r*100
impvt3 = -(error_v1_o2_r - error_SO_r)/error_SO_r*100
impvt4 = -(error_v2_r - error_SO_r)/error_SO_r*100
print('Error on {:s}, SO nor :{:4.2e}, SO {:4.2e} ({:3.1f}%), v1 opt1 {:4.2e} ({:3.1f}%), v1 opt2 {:4.2e} ({:3.1f}%), v2 {:4.2e} ({:3.1f}%) \n'.format(param,error_nor,error_SO_r,impvt1,error_v1_o1_r,impvt2,error_v1_o2_r,impvt3,error_v2_r,impvt4))










