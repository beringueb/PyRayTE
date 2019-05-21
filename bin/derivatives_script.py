import numpy as np
import matplotlib.pyplot as plt
import pickle 
lmax = 4999


with open('../data/derivatives/109A_s_nor.pick', "rb") as f:
    deriv_TT_As_nol = pickle.load(f)[:,0,0]
    
with open('../data/derivatives/tau_nor.pick', "rb") as f:
    deriv_TT_tau_nol = pickle.load(f)[:,0,0]
    
with open('../data/derivatives/109A_s.pick', "rb") as f:
    deriv_TT_As_l = pickle.load(f)[:,0,0]
    #deriv_PP_As_l = pickle.load(f)[:,-1,-1] 
    
with open('../data/derivatives/tau.pick', "rb") as f:
    deriv_TT_tau_l = pickle.load(f)[:,0,0]
    #deriv_PP_tau_l = pickle.load(f)[:,-1,-1] 
    
ell = np.linspace(2,lmax,lmax-1)
plt.figure()
plt.plot(ell, deriv_TT_As_nol)
plt.figure()
plt.plot(ell, deriv_TT_As_l)
plt.figure()
plt.plot(ell, deriv_TT_tau_nol)
plt.figure()
plt.plot(ell, deriv_TT_tau_l)

plt.show()