#!/usr/bin/python3 
""" Module that contains the definition of the fisher_matrix class and methods in it"""

import numpy as np
import os
import time
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.ticker as ticker
import pickle

class FisherMatrix():  
    """ Class that defines a fisher_matrix and methods useful to generate it, add two of them, reshuffle them, etc ... """
    
    def __init__(self,param_list,experiment):
        """ Initilaization of the fisher matrices with useful parameters
            - param_list :
            - experiment: 
        """ #TO COMPLETE
        self.experiment = experiment # experiment corresponding to the fisher matrix
        self.param_list = param_list # list of the parameters to include in the fisher matrix
        #Give a name to the fisher matrix : TT + TE + EE + PP + Rayleigh
        name = ': TT ({:d}) + '.format(experiment.lmax_T)
        if experiment.include_pol:
            name += 'TE + EE ({:d}) + '.format(experiment.lmax_P)
        if experiment.include_lens:
            name += 'PP + '
        if experiment.include_rayleigh:
            name += 'Rayleigh + {}  '.format(str(experiment.freqs))
        self.name = experiment.name + name[0:-2]
        self.fsky = experiment.fsky
        self.fisher = np.zeros(( len(self.param_list), len(self.param_list) )) # fisher matrix per se
        
    def write_to_txt(self, root, name):
        """ Method to write the fisher matrix to a txt file in root
            - name
        """
        str0 = name
        str1 = 'p' if self.experiment.include_pol else 'nop'
        str2 = 'l' if self.experiment.include_lens else 'nol'
        str3 = 'r' if self.experiment.include_rayleigh else 'nor'
        file_name = os.path.join(root,"fisher_{}_{}_{}_{}.fish".format(str0.strip(),str1,str2,str3)) 
        header = '{} \n {}'.format(self.name,str(self.param_list))
        print("Saving fisher matrix to {} ... ".format(file_name),end = '')
        np.savetxt(file_name,self.fisher,header=header)
        print("Done !")
        
        
    def get_fisher(self,setup):
        """ Method to compute the fisher matrix :
            - setup:
            -experiment:
        """
        # TO COMPLETE
        
        data_root = os.path.join(setup.data_root,self.experiment.name)
        lmax = max(self.experiment.lmax_T,self.experiment.lmax_P)
        l = np.linspace(2,lmax,lmax-1)
        #Get Cl for the fiducial cosmology
        cls_fid = read_cl(data_root, lmax, len(self.experiment.freqs), 'fiducial')
        cov_fid = compute_covariance(self.experiment,cls_fid)
        deriv = {}
        for parameter in self.param_list:
            deriv[parameter] = derivative(setup,self.experiment,data_root,parameter)
        i = 0
        time_start = time.time()
        print("Computing fisher matrix of {} for {:d} parameters ... ".format(self.name,len(self.param_list)))
        for param_i in self.param_list:
            j = 0
            for param_j in self.param_list:
                self.fisher[i,j] = self.fsky * np.sum((2*l+1.)/2.*fish_trace(lmax,deriv[param_i],deriv[param_j],cov_fid)) 
                j += 1
            i += 1
            time_tmp = time.time()
            ETA = (len(self.param_list) - i)*(time_tmp - time_start) / i
            print("{:3.1f}% done, ETA : {:2.0f} min {:2.0f} secs".format(i/len(self.param_list) * 100, ETA // 60, ETA % 60), end = "\r" )
            
        print("Done in {:2.0f} min {:2.0f} secs             ".format((time_tmp - time_start) // 60, (time_tmp - time_start) % 60))
        
    def reshuffle(self,new_param_list):
        """ Method to reshuffle the fisher matrix with a new parameter list 
            - new_param_list : New order of the parameter list, has to be of the same length than old one, separate function to fix parameter or marginalize.
        """
        old_param_list = self.param_list
        try :
            assert set(new_param_list).issubset(old_param_list)
        except AssertionError:
            print("The new parameter list needs to contain the same parameters than the old one !")
        indices = [old_param_list.index(param) for param in new_param_list]
        tmp = self.fisher
        self.fisher = tmp[np.ix_(indices, indices)]
        self.param_list = new_param_list
        return self
    
    def get_cov(self):
        """ Method to get the covariance matrix from the fisher matrix. Assumes unsused parameters have already been fixed (by reshuffling for example)"""
        cov = np.linalg.inv(self.fisher)
        return cov
        
    def __add__(self,fisher_matrix):
        """ Methos to add two fisher matrices together, reshuffle the second one to get the paramter order, change the name and add the matrices """
        try:
            assert isinstance(fisher_matrix,FisherMatrix)
        except AssertionError:
            print('Can only add two fisher matrices together')
        else:
            new_fisher = FisherMatrix(self.param_list,self.experiment)
            new_fisher.name += 'on {:2.0f}% of the sky + {} on {:2.0f}% of the sky'.format(self.fsky*100,fisher_matrix.name,fisher_matrix.fsky*100)
            fisher_matrix.reshuffle(self.param_list)
            new_fisher.fisher = self.fisher + fisher_matrix.fisher
            return new_fisher

        



def read_cl(data_root, lmax, n_freqs, parameter, direction=None):
    """ Function that gets the power spectra for a given parameter and a given direction. Uses pandas to speed up things.
        - data_root : directory containing the power pectra.
        - n_freqs : number of frequencies used by the experiment
        - parameter : parameter that we are interested in or 'fiducial' if we want to get fiducial power spectra
        - direction : left or right. If None, then read fiducial power spectra.
        returns : - cls (lmax-1 * *n_freqs + 1 * n_freqs + 1 * 5 (TT+TE+EE+TP+PP) ) TP and PP data are only taken for primary channels (index 0)
    """
    pandas = True
    try:
        import pandas as pd
    except ImportError:
        print("pandas not found, using numpy instead")
        pandas = False
    
    if direction is not None:
        filename = os.path.join(data_root,"{}_{}_".format(parameter,direction)) 
    else:
        filename = os.path.join(data_root,"{}_".format("fiducial"))
    cls = np.zeros((lmax-1,n_freqs+1,n_freqs+1,5))
    for i in range(n_freqs+1):
        for j in range(n_freqs+1):
            if pandas:
                data_read = pd.read_csv(filename + 'scalCls_lensed_{:d}_{:d}'.format(i+1,j+1), header = None, skiprows = 0, sep = '\s+').values
            else:
                data_read = np.loadtxt(filename + 'scalCls_lensed_{:d}_{:d}'.format(i+1,j+1))
            cls[:,i,j,0] = data_read[0:lmax-1,1] # TT all spectra start at ell = 2
            cls[:,i,j,1] = data_read[0:lmax-1,2] # EE
            cls[:,i,j,2] = data_read[0:lmax-1,4] # TE
    if pandas:
        data_read_phi = pd.read_csv(filename + 'scalCls_lenspot.dat', header = None, skiprows = 1, sep = '\s+').values # reading lensing potential values, files have a header
    else:
        data_read_phi = np.loadtxt(filename + 'scalCls_lenspot.dat')
    cls[:,0,0,0] = data_read_phi[0:lmax-1,1] # TT all spectra start at ell = 2
    cls[:,0,0,1] = data_read_phi[0:lmax-1,2] # EE
    cls[:,0,0,2] = data_read_phi[0:lmax-1,4] # TE
    cls[:,0,0,3] = data_read_phi[0:lmax-1,5] # PP only the primary channel is completed (for now)
    cls[:,0,0,4] = data_read_phi[0:lmax-1,6] # TP
    return cls

def kron_delta(i,j):
    """ Kronecker delta symbol, returns 1 if i==j, 0 otherwise """
    if i == j:
        return 1
    else:
        return 0
    
def compute_covariance(experiment,cls):
    """ Function to compute the covariance matrix given set of cls and an experiment
        - experiment : corresponding CMB experiment, used to know whether lensing, rayleigh or polarization should be included
        - cls : cls used to populate the covariance matrix
        return : - cov : covariance matrix
    """
    lmax = max(experiment.lmax_T,experiment.lmax_P)
    n_freqs = len(experiment.freqs)
    l_max = experiment.lmax
    # Matrix is constructed with lensing and lensing is removed at the end if unwanted
    if experiment.include_rayleigh:
        if experiment.include_pol:
            cov = np.zeros((lmax-1,2*(n_freqs+1)+1,2*(n_freqs+1)+1))
            for i in range(n_freqs+1):
                for j in range(n_freqs+1):
                    cov[:,i,j] = cls[:,0,0,0] + kron_delta(i,j)*experiment.NlTT[:,i+1] + (1.-kron_delta(i,0))*cls[:,i,0,0] + (1.-kron_delta(j,0))*cls[:,0,j,0] + (1-kron_delta(i,0))*(1-kron_delta(j,0))*cls[:,i,j,0] # TT part
                    cov[:,i+n_freqs+1,j+n_freqs+1] = cls[:,0,0,1] + kron_delta(i,j)*experiment.NlEE[:,i+1] + (1.-kron_delta(i,0))*cls[:,i,0,1] + (1.-kron_delta(j,0))*cls[:,0,j,1] + (1-kron_delta(i,0))*(1-kron_delta(j,0))*cls[:,i,j,1] # EE part
                    cov[:,i+n_freqs+1,j] = cls[:,0,0,2] + (1.-kron_delta(i,0))*cls[:,i,0,2] + (1.-kron_delta(j,0))*cls[:,0,j,2] + (1-kron_delta(i,0))*(1-kron_delta(j,0))*cls[:,i,j,2] # TE part
                    cov[:,i,j+n_freqs+1] = cls[:,0,0,2] + (1.-kron_delta(i,0))*cls[:,i,0,2] + (1.-kron_delta(j,0))*cls[:,0,j,2] + (1-kron_delta(i,0))*(1-kron_delta(j,0))*cls[:,i,j,2] # TE part
        else:
            cov = np.zeros((l_max-1,(n_freqs+1)+1,(n_freqs+1)+1))
            for i in range(n_freqs+1):
                for j in range(n_freqs +1):
                    cov[:,i,j] = cls[:,0,0,0] + kron_delta(i,j)*experiment.NlTT[:,i+1] + (1.-kron_delta(i,0))*(cls[:,i,0,0]) + (1.-kron_delta(j,0))*(cls[:,0,j,0]) + (1-kron_delta(i,0))*(1-kron_delta(j,0))*(cls[:,i,j,0]) # TT part
    else:
        if experiment.include_pol:
            cov = np.zeros((lmax-1,3,3))
            cov[:,0,0] = cls[:,0,0,0] + experiment.NlTT[:,1] #+ 3.
            cov[:,0,1] = cls[:,0,0,2]
            cov[:,1,0] = cls[:,0,0,2]
            cov[:,1,1] = cls[:,0,0,1] + experiment.NlEE[:,1]
        else:
            cov = np.zeros((lmax-1,2,2))
            cov[:,0,0] = cls[:,0,0,0] + experiment.NlTT[:,1]
    if experiment.include_lens:
        cov[:,-1,0] = 0*cls[:,0,0,4] # TP part
        cov[:,0,-1] = 0*cls[:,0,0,4] # TP part
        cov[:,-1,-1] = cls[:,0,0,3] + experiment.NlPP[:,1] #PP part
        return cov
    else:
        cov_out = np.delete(cov,(-1),axis = 1)
        cov_out_out = np.delete(cov_out,(-1),axis = 2)
        return cov_out_out
        

def derivative(setup,experiment,data_root,parameter):
    """ Function to compute the derivative of the the covariance matrix with respect to a given parameter.
        - setup :  
        - experiment : corresponding CMB experiment, used to get the number of frequency channels, lmax
        - data_root : directory containing the powers pectra
        - parameter : parameter with respect to which take the derivative
        return : deriv
    """ #TO COMPLETE
    lmax = max(experiment.lmax_T,experiment.lmax_P)
    n_freqs = len(experiment.freqs)
    #Reading data in every direction
    cls_right = read_cl(data_root, lmax, n_freqs, parameter, direction='right')
    cls_left = read_cl(data_root, lmax, n_freqs, parameter, direction='left')
    #Compute the corresponding covariance matrices
    cov_right = compute_covariance(experiment,cls_right)
    cov_left = compute_covariance(experiment,cls_left)
    #Compute the derivative (2pts for now, maybe update at some point ?)
    deriv = (cov_right - cov_left)/(2.*setup.step[parameter])
    file_name = '../data/derivatives/{:s}_nol_r.pick'.format(parameter)
    with open(file_name, 'wb') as f:
        pickle.dump(deriv,f)
    return deriv
    
def fish_trace(lmax,deriv_i,deriv_j,cov_fid):
    """ Function that comput the trace part of the fisher formula.
        - lmax : ell max
        - deriv_i : derivative of the covariance matrix with respect to parameter_i
        - deriv_j : derivative of the covariance matrix with respect to parameter_j
        - cov_fid : fiducial covariance matrix
        return : tr(inv(cov_fid)*deriv_i*inv(cov_fid)*deriv_j) array of lmax-1 entries
    """
    
    tr = np.zeros(lmax-1)
    for ell in range(lmax-1):
        inv = np.linalg.inv(cov_fid[ell,:,:])
        tr[ell] = np.trace(np.dot(inv,np.dot(deriv_i[ell,:,:],np.dot(inv,deriv_j[ell,:,:]))))
    return tr
    
def plot_cov_ellipse(cov, pos, nstd=2, ax=None, **kwargs):
    if ax is None:
        ax = plt.gca()
    a = np.sqrt((cov[0,0]+cov[1,1])/2. + np.sqrt((cov[0,0]-cov[1,1])**2/4. + cov[0,1]**2))
    b = np.sqrt((cov[0,0]+cov[1,1])/2. - np.sqrt((cov[0,0]-cov[1,1])**2/4. + cov[0,1]**2))
    theta = np.degrees(np.arctan2(2*cov[0,1],cov[0,0]-cov[1,1]))/2
    alphas = np.asarray((2.3,6.17,11.8))
    alpha = np.sqrt(alphas[nstd-1]) 
    # Width and height are "full" widths, not radius
    ellip = Ellipse(xy=pos, width=2*alpha*a, height=2*alpha*b, angle=theta, **kwargs)
    ax.add_artist(ellip)
    return ellip
    
  
def print_errors(list_fisher):
    """ Print marginalized errors on paramters. Assumed fisher matrices in list_fisher are reshuffled to have the same parameter order.
        - list_fisher : list of the fisher matrices from which to get the errors.
        return : string
    """
    cov_list = [fish.get_cov() for fish in list_fisher]
    parameter_list = list_fisher[0].param_list
    string = ''
    j = 0
    for param in parameter_list:
        string += " Error on {}   :  ".format(param)
        for i in range(len(cov_list)):
            string += '{:5.3e}   |'.format(np.sqrt(cov_list[i][j,j]) )
        string += '\n'
        j +=1
    return string
    
def get_error(fisher,parameter):
    """ Get the error on a given parameter marginalized over every other parameter in the fisher matrix
        - fisher : fisher matrix
        - parameter : parameter we are interested in
        return : error : error on parameter marginalized over every other paramters in the fisher matrix
    """
    try : 
        index = fisher.param_list.index(parameter) 
    except ValueError:
        print('{} is not in the list of parameters for the fiser matrix !'.format(parameter))
        error = np.inf
    else:
        error = np.sqrt(np.linalg.inv(fisher.fisher)[index,index])
    return error
        
        

def plot_triangular(setup, list_fisher, parameter_list=None):
    """ Plot triangular plot for list of fisher matrix
        - setup_used, must be the same for every fisher matrices, at least in terms of fiducial values and step sizes
        - list_fisher : list of fisher matrices to plot
        - parameter_list : list of parameters to plot, if None uses the paramter list of the first fisher matrix in the list
    """


    if parameter_list is None:
        parameter_list = list_fisher[0].param_list
    list_to_plot = []
    for fish in list_fisher:
        list_to_plot.append(fish.reshuffle(parameter_list))
    n_params = len(parameter_list)
    color = ('k','c','b','y','g','r','gold','sienna','coral','navy')
    names = {'H0':r'$H_0$','A_s':r'$A_s$','109A_s':r'$10^9 A_s$','ln10A_s':r'$\ln(10^{10}a_s)$','N_eff':r'$N_{eff}$','Y_He':r'$Y_{He}$','ombh2':r'$\Omega_bh^2$','omch2':r'$\Omega_ch^2$','theta_MC':r'$\theta_{MC}$','n_s':r'$n_s$','tau':r'$\tau$'}
    center = [setup.fiducial[param] for param in parameter_list]
    delta = [setup.step[param] for param in parameter_list]
    delta[2] /= 8
    delta[0] /= 3
    label_list = ['SO + PLANCK : TT + TE + EE + PP + Rayleigh', 'SO + PLANCK : TT + TE + EE + PP', 'PLANCK : TT + TE + EE + PP' ]
    plt.rc('text', usetex = True) 
    fig,axs = plt.subplots(n_params,n_params)
    plt.subplots_adjust(wspace=0, hspace=0,left=0.1,right=0.98,top = 0.98,bottom=0.1)
    for i in range(n_params):
        for j in range(i+1):
            if i != j:
                for k in range(len(list_to_plot)):
                    cov = list_to_plot[k].get_cov()[np.ix_([j,i],[j,i])]
                    e = plot_cov_ellipse(cov,[center[j],center[i]],ax=axs[i,j],nstd = 1,lw = 1.2, ec = color[k], fc = 'None',alpha = 1,label = label_list[k])
                
                axs[i,j].scatter(center[j],center[i], c = 'r', marker = '+')
                axs[i,j].grid(color='k', linestyle = '-.', linewidth = 0.5, alpha = 0.5)
                axs[i,j].set_xlim(center[j]-1*delta[j], center[j]+1*delta[j])
                axs[i,j].set_ylim(center[i]-1*delta[i], center[i]+1*delta[i])
                axs[i,j].xaxis.set_ticks(np.linspace(center[j]-1*delta[j], center[j]+1*delta[j],9))
                axs[i,j].yaxis.set_ticks(np.linspace(center[i]-1*delta[i], center[i]+1*delta[i],7))
                for label in axs[i,j].xaxis.get_ticklabels()[::2]:
                    label.set_visible(False)
                for label in axs[i,j].yaxis.get_ticklabels()[::2]:
                    label.set_visible(False)
                fig.delaxes(axs[j,i]) #delete unused plots
                
            else: 
                xx = np.linspace(center[i]-6*delta[i], center[i]+6*delta[i], 2000)
                for k in range(len(list_to_plot)):
                    cov = list_to_plot[k].get_cov()
                    yy = np.exp(-(xx-center[i])**2/(2*cov[i,i]))
                    axs[i,j].plot(xx,yy,c = color[k], lw = 1.2, label = label_list[k])
                axs[i,j].grid(color='k', linestyle = '-.', linewidth = 0.5, alpha = 0.5)
                axs[i,j].xaxis.set_ticks(np.linspace(center[j]-1*delta[j],center[j]+1*delta[j],9))
                axs[i,j].set_xlim(center[j]-1*delta[j], center[j]+1*delta[j])
                axs[i,j].set_ylim(0,1)
                for label in axs[i,j].xaxis.get_ticklabels()[::2]:
                    label.set_visible(False)
            axs[i,j].tick_params('both', labelsize = 7)
            axs[0,0].set_yticklabels([])
            axs[n_params-1,j].set_xlabel(names[parameter_list[j]], fontsize = 16)
            x_formatter = ticker.ScalarFormatter(useOffset=True)
            x_formatter.set_scientific(True)
            axs[n_params-1,j].xaxis.set_major_formatter(x_formatter)
            axs[0,0].legend(loc = 'center',prop = {'size':16},fancybox = True, shadow = True,bbox_to_anchor = (3.5,0.45))
            if i != n_params-1 :
                axs[i,j].set_xticklabels([])
            if j != 0 :
                axs[i,j].set_yticklabels([])
            if i != 0:
                axs[i,0].set_ylabel(names[parameter_list[i]], fontsize = 16)
                y_formatter = ticker.ScalarFormatter(useOffset=True)
                y_formatter.set_scientific(True)
                axs[i,0].yaxis.set_major_formatter(y_formatter)
    plt.show()
                

                

                    
                    
    
                         
                         
                         
                         
                        
    
    
    
    
        
    
        
    
           
    
             
   
   
        
        
        
        
        
        
        
        
        
        

