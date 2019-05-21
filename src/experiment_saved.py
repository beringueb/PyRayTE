#!/usr/bin/python3

"""Modules that contains classes and methods to read a .ini file (in ../ini/) and initialize several classes. """

import os
import configparser
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


class Experiment() : 
    """ Class that contains information about a given CMB experiment : 
        - name : name of the experiment
        - include_pol : whether to include polarization data from this experiment
        - include_lens : whether to include lensing potential information from this experiment
        - include_rayleigh : whether to include to Rayleigh scattering signal
        - freqs : list of frequency bands for this experiment (does not include primary channel)
        - noise_file : location of the noise power spectra for this experiment, template is given in "../data/noise/template.dat"
        - fsky : fraction of the sky covered by the experiment
        - lmax_T : ell up to which include temperature data
        - lmax_P : ell up to which include polarization data
        - lmin : ell from which to include data (assumed to the same in temperature and polarization)
        - NlTT : (lmax-1,n_freqs+2) temperature noise power spectra. Columns are ell, primary channel (combined freqs), freqs bands, ell from 2 to lmax, units are ell*(ell+1.)/2pi*N_l
        - NlEE : (lmax-1,n_freqs+2) polarization noise power spectra. Columns are ell, primary channel (combined freqs), freqs bands, ell from 2 to lmax, units are ell*(ell+1.)/2pi*N_l
        - NlPP : (lmax-1,2) lensing noise power spectra. Columns are ell, primary channel (combined freqs), ell from 0 to lmax, units are (ell*(ell+1.))**2/2pi*N_l
    """
    def __init__(self,name,include_pol,include_lens,include_rayleigh,freqs,noise_file,fsky,lmax_T,lmax_P,lmin):
        """Initialization of an experiment"""
        lmax = max(lmax_T,lmax_P)
        self.name = name
        self.include_pol = include_pol
        self.include_lens = include_lens
        self.include_rayleigh = include_rayleigh
        self.freqs = freqs
        self.fsky = fsky
        self.lmax_T = lmax_T
        self.lmax_P = lmax_P
        self.lmax = lmax
        self.lmin = lmin
        self.NlTT = np.zeros((lmax-1,len(freqs)+2)) + 1E15
        self.NlEE = np.zeros((lmax-1,len(freqs)+2)) + 1E15
        self.NlPP = np.zeros((lmax-1,2)) + 1E15
        ell = np.linspace(2,lmax,lmax-1)
        self.NlTT[:,0] = ell
        self.NlEE[:,0] = ell
        self.NlPP[:,0] = ell
        self.noise_file = noise_file
            
        
    def read_noise(self):
        """ Reading noise power spectra. Assumes file already exists."""
        print("Reading noise power spectra from {} ... ".format(self.noise_file), end = "")
        noise_read = pd.read_csv(self.noise_file, sep='\s+', skiprows = [0,1,2], header = None )
        noise = noise_read.values 
        self.NlTT[:,2:] = noise[2:self.lmax+1,1:(len(self.freqs)+1)]  ## MODIFIED 2: !!!!!!!
        self.NlEE[:,2:] = noise[2:self.lmax+1,(len(self.freqs)+1):-1]
        if (np.shape(noise)[1]-1) % 2 != 0 : 
            self.NlPP[:,1] = noise[0:self.lmax-1,-1] # in that case, no lensing noise is available for this experiment so leave it to inf
        print("Done ! ")

    def max_min(self):
        """Set noise to infinity for ell outside range"""
        
        #Deleting contributions from ell < lmin
        
        self.NlTT[0:self.lmin-2,1:] = 1E15 
        self.NlEE[0:self.lmin-2,1:] = 1E15
        self.NlPP[0:self.lmin-2,1] = 1E15
        
        #Deleting contributions from ell > lmax
        
        self.NlTT[self.lmax_T:, 1:] = 1E15
        self.NlEE[self.lmax_P:, 1:] = 1E15
        if self.include_lens:
            self.NlPP[min(self.lmax_P,self.lmax_T):, 1] = 1E15
        else:
            self.NlPP[:,1] = 1E15
        if not self.include_pol:
            self.NlEE[:,1:] = 1E15

        

    def combine_primary(self, list_freqs = 'all'):
        """ Combine frequency channels to get the the primary CMB noise
            if list_freqs is all then all frequency bands are used
            otherwise uses frequencies in list_freqs
        """

        try : 
            assert isinstance(list_freqs,list) or list_freqs == 'all'
        except AssertionError : 
            print("list_freqs must be either a list (not a tuple) or 'all' ... ")
        lmax = self.lmax
        temp_T = np.zeros(lmax-1) + 1E-25
        temp_E = np.zeros(lmax-1) + 1E-25
        if list_freqs == 'all':
            list_freqs = self.freqs
        print("Combining {} GHz channels into primary CMB noise ... ".format(str(list_freqs)), end = "")
        for freq in list_freqs: #simple inverse variance weighting of the different frequency channels, could think of the possibility to allow different weights
            index = self.freqs.index(freq)
            temp_T += 1./self.NlTT[:,index+2]
            temp_E += 1./self.NlEE[:,index+2]

        self.NlTT[:,1] = 1./temp_T
        self.NlEE[:,1] = 1./temp_E
        print("Done !")
        
    def copy(self):
        """ Method to return a copy of on an experiment"""
        copy = Experiment(self.name,self.include_pol,self.include_lens,self.include_rayleigh,self.freqs,self.noise_file,self.fsky,self.lmax_T,self.lmax_P,self.lmin)
        copy.lmax = self.lmax 
        copy.NlTT = self.NlTT 
        copy.NlEE = self.NlEE 
        copy.NlPP = self.NlPP
        return copy
        
        
    def plot(self,list_freqs = None, factor = True):
        """ Method to plot TT and EE noise curves.
            - list_freqs : if None plot primary channels, else plot frequencies from list_freqs + primary
            - factor : if True, plot ell*(ell+1)/2pi * N_l
        """
        import matplotlib.pyplot as plt
        ell = self.NlTT[:,0]
        if factor :
            alpha = 1
        else :
            alpha = 2*np.pi/ell/(ell+1)
        f, axs = plt.subplots(1,2, sharex = True, sharey = True)
        color = ['k','c','b','y','g','r','gold','sienna','coral','navy','darkblue','brown','lime','maroon','olive','salmon']
        if list_freqs is None:
            axs[0].loglog(ell,self.NlTT[:,1]*alpha, label = "Primary channel")
            axs[1].loglog(ell,self.NlEE[:,1]*alpha, label = "Primary channel")
        else :
            k = 0    
            for i in range(len(list_freqs)+1):
                fr = list_freqs[i]
                index = self.freqs.index(fr)
                axs[0].loglog(ell,self.NlTT[:,2+index]*alpha, c = color[k], label = "{:d} GHz".format(fr))
                axs[1].loglog(ell,self.NlEE[:,2+index]*alpha, c = color[k], label = "{:d} GHz".format(fr))
                k+=1
        axs[0].set_title("N($\ell$) Temperature")
        axs[1].set_title("N($\ell$) Polarization")
        axs[0].set_ylabel("N($\ell$) [$\mu K^2$-SR] " )
        axs[1].set_ylabel("N($\ell$) [$\mu K^2$-SR] " )
        axs[0].set_xlabel("$\ell$")
        axs[1].set_xlabel("$\ell$")
        axs[0].set_ylim(1e-8,10)
        axs[0].legend(bbox_to_anchor=(0.02, 0.80,0.96,0.90), loc=3,ncol=4, mode="expand", borderaxespad=0.)
        axs[1].legend(bbox_to_anchor=(0.02, 0.80,0.96,0.90), loc=3,ncol=4, mode="expand", borderaxespad=0.)
        axs[0].set_xlim(10, int(ell[-1]))
        plt.suptitle("Noise power spectra for {}".format(self.name))
        plt.show()
        
def combine(list_experiments):
    """ Method to combine experiments. On same patch of the sky, inverse variance weight primary channels and add both lists of frequencies.
        Therefore it returns a list of combined experiments on independent patches of the sky. """
    fsky_list = list(set([experiment.fsky for experiment in list_experiments]))
    fsky_list.sort()
    N_patches = len(fsky_list) # number of independent experiments (on non overlapping part of the sky)
    experiments_to_include = list_experiments.copy()
    list_patches =[]
    fsky_current = 0.
    for i in range(N_patches):
        fsky = fsky_list[i]
        list_patches.append([expe.copy() for expe in experiments_to_include])
        for experiment in list_patches[i] :
            experiment.fsky = fsky
        for experiment in experiments_to_include:
            if experiment.fsky-fsky_current == fsky:
                experiments_to_include.remove(experiment)
        fsky_list = [i - fsky for i in fsky_list]
        fsky_current = fsky
    combined_experiments = []
    for patch in list_patches:
        combined_experiments.append(combine_on_same_patch(patch))
    return combined_experiments
            
        
def combine_on_same_patch(list_experiments):
    """ Function that combines experiments on the same patch of the sky. Noise spectra are inverse variance weighted."""
    N_experiments = len(list_experiments)
    if N_experiments == 1:
        print('{} on {:4.1f}% of the sky'.format(list_experiments[0].name,list_experiments[0].fsky*100))
        #print(list_experiments[0].NlTT,list_experiments[0].NlEE,list_experiments[0].NlPP, sep = '\n')
        return list_experiments[0]
    name = list_experiments[0].name
    for i in range(N_experiments-1):
        name += ' + {}'.format(list_experiments[i+1].name)
    include_pol = any([expe.include_pol for expe in list_experiments])
    include_lens = any([expe.include_lens for expe in list_experiments])
    include_rayleigh= any([expe.include_rayleigh for expe in list_experiments])
    fsky = list_experiments[0].fsky #at this point all experiments in the list should have the same fsky
    freqs = []
    for expe in list_experiments:
        if expe.include_rayleigh:
            freqs += expe.freqs
    freqs.sort()
    lmax_T = max([expe.lmax_T for expe in list_experiments])
    lmax_P = max([expe.lmax_P for expe in list_experiments])
    lmin = max([expe.lmin for expe in list_experiments])
    
    combined = Experiment(name.replace(' ',''),include_pol,include_lens,include_rayleigh,freqs,'',fsky,lmax_T,lmax_P,lmin)
    lmax = combined.lmax

    NlTT = np.zeros((lmax-1,2 + len(combined.freqs))) 
    NlEE = np.zeros((lmax-1,2 + len(combined.freqs))) 
    NlPP = np.zeros((lmax-1,2)) 
    i = 0
    for freq in combined.freqs:
        noiseTT_tmp = np.zeros(lmax-1) + 1e-15
        noiseEE_tmp = np.zeros(lmax-1) + 1e-15
        #look in each experiment if there is a channel at that frequency 
        for expe in list_experiments:
            try:
                ind = expe.freqs.index(freq)
            except ValueError:
                pass #freq not in this experiment ... 
            else:
                noiseTT_tmp[0:expe.lmax-1] += 1./(expe.NlTT[:,2+ind]) #scaling with fsky self.fsky/expe.fsky
                noiseEE_tmp[0:expe.lmax-1] += 1./(expe.NlEE[:,2+ind])
        NlTT[:,2+i] = 1./noiseTT_tmp
        NlEE[:,2+i] = 1./noiseEE_tmp
        i+=1
    noiseTT_tmp = np.zeros(lmax-1) + 1E-15
    noiseEE_tmp = np.zeros(lmax-1) + 1E-15
    noisePP_tmp = np.zeros(lmax-1) + 1E-15    
    for expe in list_experiments: #same for primary channels
        noiseTT_tmp[0:expe.lmax-1] += 1./(expe.NlTT[:,1])
        noiseEE_tmp[0:expe.lmax-1] += 1./(expe.NlEE[:,1])
        noisePP_tmp[0:expe.lmax-1] += 1./(expe.NlPP[:,1])
    NlTT[:,1] = 1./noiseTT_tmp
    NlEE[:,1] = 1./noiseEE_tmp
    NlPP[:,1] = 1./noisePP_tmp
    NlTT[:,0] = np.linspace(2,lmax,lmax-1)
    NlEE[:,0] = np.linspace(2,lmax,lmax-1)
    NlPP[:,0] = np.linspace(2,lmax,lmax-1)
    combined.NlTT = NlTT
    combined.NlEE = NlEE
    combined.NlPP = NlPP
    print('{} on {:4.1f}% of the sky'.format(combined.name,combined.fsky*100))
    return combined
    
    

class Setup():
    """ Class that summarizes the setup used for the Fisher Forecast
        - param_list : lsit of paramters used for Fisher analysis
        - data_root : directory where power spectra are written
        - fiducial : dict with the fiducial values for the parameters used
        - step : dict of step sizes taken for the derivatives
        - list_experiments : list of (combined) experiments used in the analysis
        - use_BBN : 1 if we want to use BBN consistency for Y_He. If fixed to a given value, assume this value for y_He. Overiden if y_He is in param_list
        - mass_neutrinos : mass of the neutinos (in eV)
        - lmax : ell max used in the analysis, used by CAMB but doesn't override lmax for experiments
        - l_boost : CAMB parameter, keep more terms in the hierarchy evolution
        - boost : CAMB parameter, increase accuracy_boost to decrease time steps, use more k values
    """
    def __init__(self,param_list,data_root, use_BBN, mass_neutrinos = 60., l_boost = 1, boost = 1):
        """ Initialization """
        self.param_list = param_list
        self.data_root = data_root
        self.use_BBN = use_BBN
        self.mass_neutrinos = mass_neutrinos
        self.l_boost = l_boost
        self.boost = boost
        
    def get_list_experiments(self, list_experiments):
        """ Method to get the list of combined experiment in the setup. Separated from the rest of the __init__ method to make parsing easier"""
        self.list_experiments = list_experiments
        self.lmax = max([expe.lmax for expe in list_experiments])
    
    def get_fiducial(self, fid_file):
        """ Read fiducial values and step sizes for parameters from an external file 
            - fid_file : file from which values are read. Lines must be : name | fiducial value | step size with one header line starting with #
        """  
        fid = {}
        step = {}
        with open(fid_file,'r') as f:
            for line in f:
                 if '####' in line:
                    break
                 else:  
                    name,fiducial_value,step_size = line.split('|')
                    fid[name.strip()] = float(fiducial_value)
                    step[name.strip()] = float(step_size)
        try:
            assert set(self.param_list).issubset(fid.keys())
        except AssertionError:
            print("Couldn't load fiducial values and step sizes for some parameters. Check {} !".format(fid_file))
        
        self.fiducial = fid
        self.step = step
        
def parser(inifile):
    """ Function that parses arguments from .ini file
        - inifile : loaction of the .ini file
    """
    config = configparser.SafeConfigParser()
    config.read(inifile)        
    sections = config.sections()    
    
    ## READING SETUP CONFIGURATION ##
    setup_sec = sections[0]
    param_list = [param for param in config.get(setup_sec,'parameters').split(',')]
    print('Fisher analysis will run for {}'.format(param_list))
    try:
        assert set(param_list).issubset(['H0','A_s','109A_s','ln10A_s','N_eff','Y_He','ombh2','omch2','theta_MC','n_s','tau','mass_neutrinos'])
    except AssertionError:
        print("Some parameters are invalid. Check you ini file")
    data_root = config.get(setup_sec,'data_root')
    use_BBN = config.getfloat(setup_sec,'use_BBN')
    if config.has_option(setup_sec,'mass_neutrinos') :
        mass_neutrinos = config.getfloat(setup_sec,'mass_neutrinos')
    else:
        mass_neutrinos = 60.
    if config.has_option(setup_sec,'CAMB_l_boost'):
        l_boost = config.getfloat(setup_sec,'CAMB_l_boost')
    else:
        l_boost = 1
    if config.has_option(setup_sec,'CAMB_accuracy_boost'):
        boost = config.getfloat(setup_sec,'CAMB_accuracy_boost')
    else:
        boost = 1    
    
    setup = Setup(param_list,data_root, use_BBN, mass_neutrinos, l_boost, boost)
    list_experiments = []
       
    ## READING EXPERIMENTS CONFIGURATIONS ##
    for experiment in sections[1:]:
        name = config.get(experiment,'name')
        fsky = config.getfloat(experiment,'fsky')
        include = config.getboolean(experiment,'include')
        include_pol = config.getboolean(experiment,'polarization')
        include_lens = config.getboolean(experiment,'lensing')
        include_rayleigh = config.getboolean(experiment,'rayleigh')
        lmin = config.getint(experiment,'lmin')
        lmax_T = config.getint(experiment,'lmax_T')
        lmax_P = config.getint(experiment,'lmax_P')
        freqs = [ int(fr) for fr in config.get(experiment,'freqs').split(',')]      
        noise_file = config.get(experiment,'noise_file')
        print(noise_file)

        expe = Experiment(name,include_pol,include_lens,include_rayleigh,freqs,noise_file,fsky,lmax_T,lmax_P,lmin)
        if include:
            list_experiments.append(expe)
        
    return setup,list_experiments
    
def ini_driver(inifile):
    """ Function that calls parser, read noise, combine experiments and returns setup"""
    setup, list_experiments = parser(inifile)
    
    for experiment in list_experiments:
        experiment.read_noise()
        experiment.max_min()
        if experiment.name == 'PICO':
            experiment.combine_primary( experiment.freqs[:3])
        elif experiment.name == 'CMB-HD':#
            experiment.combine_primary(experiment.freqs[2:3])
            #print(experiment.freqs[2:3])
        else: experiment.combine_primary()


    combined_experiments = combine(list_experiments)
    setup.get_list_experiments(combined_experiments)

    return setup
    
    
if __name__ == "__main__":
    exp1 = Experiment('pla',True,True,True,[0,1,2],'../data/noise_spectra/template.dat',0.5,100,100,2)
    exp2 = Experiment('pla2',True,True,True,[0,1,2],'../data/noise_spectra/template.dat',0.75,100,100,2)
    exp3 = Experiment('pla2',True,True,True,[0,1,2],'../data/noise_spectra/template.dat',0.80,100,100,2)
    combined = combine([exp1,exp2,exp3])
    for expe in combined:
        print(expe.name,expe.fsky, sep = '\n')

    

