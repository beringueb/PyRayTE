#!/usr/bin/python3

"""Modules that contains classes and methods to read a .ini file (in ../ini/) and initialize several classes. """

import os
import configparser
import numpy as np
import pandas as pd


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
        self.NlTT = np.zeros((lmax-1,len(freqs)+2)) + np.inf
        self.NlEE = np.zeros((lmax-1,len(freqs)+2)) + np.inf
        self.NlPP = np.zeros((lmax-1,2)) + np.inf
        ell = np.linspace(2,lmax,lmax-1)
        self.NlTT[:,0] = ell
        self.NlEE[:,0] = ell
        self.NlPP[:,0] = ell
        self.noise_file = noise_file
            
        
    def read_noise(self):
        """ Reading noise power spectra. Assumes file already exists."""
        print("Reading noise power spectra from {} ... ".format(self.noise_file), end = "")
        noise_read = pd.read_csv(self.noise_file, sep='\s+', skiprows = 3)
        noise = noise_read.values
        self.NlTT[:,2:] = noise[0:self.lmax-1,2:(len(self.freqs)+2)]
        self.NlEE[:,2:] = noise[self.lmax-1,(len(self.freqs)+2):-1]
        if (np.shape(noise)[1]-1) % 2 != 0 : 
            self.NlPP[:,:] = noise[0:self.lmax-1,-1] # in that case, no lensing noise is available for this experiment so leave it to inf
        #print(self.NlTT)
        #print(self.NlEE)
        #print(self.NlPP)
        print("Done ! ")

    def max_min(self):
        """Set noise to infinity for ell outside range"""
        
        #Deleting contributions from ell < lmin
        
        self.NlTT[0:self.lmin-2,1:] = np.inf
        self.NlEE[0:self.lmin-2,1:] = np.inf
        self.NlPP[0:self.lmin-2,1] = np.inf
        
        #Deleting contributions from ell > lmax
        
        self.NlTT[self.lmax_T:, 1:] = np.inf
        self.NlEE[self.lmax_P:, 1:] = np.inf
        if self.include_lens:
            self.NlPP[min(self.lmax_P.self.lmax_T):, 1] = np.inf
        else:
            self.NlPP[:,1] = np.inf
        

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
        temp_T = np.zeros(lmax-1)
        temp_E = np.zeros(lmax-1)
        if list_freqs == 'all':
            list_freqs = self.freqs
        print("Combining {} GHz channels into primary CMB noise ... ".format(str(list_freqs)), end = "")
        for freq in list_freqs: #simple inverse variance weighting of the different frequency channels, could think of the possibility to allow different weights
            index = list_freqs.index(freq)
            temp_T += 1./self.NlTT[:,index+2]
            temp_E += 1./self.NlEE[:,index+2]
        self.NlTT[:,1] = 1./temp_T
        self.NlEE[:,1] = 1./temp_E
        print("Done !")
        
        
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
            for fr in list_freqs:
                i = self.freqs.index(fr)
                axs[0].loglog(ell,self.NlTT[:,2+i]*alpha, c = color[k], label = "{:d} GHz".format(fr))
                axs[1].loglog(ell,self.NlEE[:,2+i]*alpha, c = color[k], label = "{:d} GHz".format(fr))
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
                
 
        
class CombinedExperiment():
    """ Class that contains conbined data for experiments observing the same sky fracion
        - list_experiments : list of the experiments to combine
        - name : name given to the combination of these experiments ('name1' + 'name2' + ...)
        - include_P : whether to include polarization data from these experiments (include_P1 or include_P2 or ...)
        - include_lens : whether to include lensing potential information from these experiments (include_lens1 or include_lens2 or ...)
        - include_rayleigh : whether to include to Rayleigh scattering signal from these experiments (include_rayleigh1 or include_rayleigh2 or ...)
        - freqs : list of frequency bands for these combined experiments, concatenation of frequency channels for experiments that include rayleigh
        - fsky : fraction of the sky covered by the experiments, min of sky fractions of considered experiments
        - lmax_T : ell up to which include temperature data, max of lmax_T of considered experiments
        - lmax_P : ell up to which include polarization data, max of lmax_P of considered experiments
        - lmin : ell from which to include data, max of lmin of considered experiments
        - NlTT : (lmax+1,n_freqs+2) temperature noise power spectra. Columns are ell, primary channel (combined freqs), freqs bands, ell from 0 to max lmax, units are ell*(ell+1.)/2pi. Inverse variance sum of NlTT for combined experiments.
        - NlEE : (lmax+1,n_freqs+2) polarization noise power spectra. Columns are ell, primary channel (combined freqs), freqs bands, ell from 0 to max lmax, units are ell*(ell+1.)/2pi.
Inverse variance sum of NlEE for combined experiments that include polarization.
        - NlPP : (lmax+1,2) lensing noise power spectra. Columns are ell, primary channel (combined freqs), ell from 0 to max lmax, units are (ell*(ell+1.))**2/2pi. Inverse variance sum of NlPP for combined experiments that include lensing.
    """
    def __init__(self, list_experiments,fsky):
        """ Initialization of a combined experimentnumpy.core._internal.AxisError: axis 1 is out of bounds for array of dimension 1
            - list_experiments : list of the experiments to combine.
            - fsky : fraction of the sky covered by all the experiments in list_experiments
        """
        self.list_experiments = list_experiments
        self.fsky = fsky
        name = ''
        for expe in list_experiments:
            name += '{} + '.format(expe.name)
        self.name = name[0:-3].replace(' ','') # name of the experiment
        print("Combining {} on {:3.1f}% of the sky.".format(self.name,self.fsky*100))
        
        include_P = False
        include_lens = False
        include_rayleigh = False
        for expe in list_experiments:
            include_P = include_P or expe.include_P
            include_lens = include_lens or expe.include_lens
            include_rayleigh = include_rayleigh or expe.include_rayleigh
        self.include_P = include_P
        self.include_lens = include_lens
        self.include_rayleigh = include_rayleigh
        
        freqs = []
        for expe in list_experiments:
            if expe.include_rayleigh:
                freqs += expe.freqs  # include the freqs of an experiment only if include Rayleigh from that experiment
        
        self.freqs = freqs
        self.freqs.sort()
        self.lmax_T = max([expe.lmax_T for expe in list_experiments])
        self.lmax_P = max([expe.lmax_P for expe in list_experiments])
        self.lmin = min([expe.lmin for expe in list_experiments])

        lmax = max(self.lmax_T,self.lmax_P)

        NlTT = np.zeros((lmax-1,2 + len(self.freqs))) 
        NlEE = np.zeros((lmax-1,2 + len(self.freqs))) 
        NlPP = np.zeros((lmax-1,2)) + np.inf
        i = 0
        for freq in self.freqs:
            noiseTT_tmp = np.zeros(lmax-1) + 1e-15
            noiseEE_tmp = np.zeros(lmax-1) + 1e-15
        
            #look in each experiment if there is a channel at that frequency 
            for expe in list_experiments:
                try:
                    ind = expe.freqs.index(freq)
                except ValueError:
                    pass #freq not in this experiment ... 
                else:
                    noiseTT_tmp[0:expe.lmax-1] += 1./(self.fsky/expe.fsky * expe.NlTT[:,2+ind]) #scaling with fsky
                    noiseEE_tmp[0:expe.lmax-1] += 1./(self.fsky/expe.fsky * expe.NlEE[:,2+ind])
            NlTT[:,2+i] = 1./noiseTT_tmp
            NlEE[:,2+i] = 1./noiseEE_tmp
            i+=1
        noiseTT_tmp = np.zeros(lmax-1) + 1E-15
        noiseEE_tmp = np.zeros(lmax-1) + 1E-15
        noisePP_tmp = np.zeros(lmax-1) + 1E-15    
        for expe in list_experiments:
            noiseTT_tmp[0:expe.lmax-1] += 1./(self.fsky/expe.fsky * expe.NlTT[:,1])
            noiseEE_tmp[0:expe.lmax-1] += 1./(self.fsky/expe.fsky * expe.NlEE[:,1])
            noisePP_tmp[0:expe.lmax-1] += 1./(self.fsky/expe.fsky * expe.NlPP[:,1])
        NlTT[:,1] = 1./noiseTT_tmp
        NlEE[:,1] = 1./noiseEE_tmp
        NlPP[:,1] = 1./noisePP_tmp
        NlTT[:,0] = np.linspace(2,lmax,lmax-1)
        NlEE[:,0] = np.linspace(2,lmax,lmax-1)
        NlPP[:,0] = np.linspace(2,lmax,lmax-1)
        self.NlTT = NlTT
        self.NlEE = NlEE
        self.NlPP = NlPP
        
        
    def plot(self,list_freqs = None):
        """ Method to plot TT and EE noise curves.
            - list_freqs : if None plot primary channels, else plot frequencues from list_freqs
        """
        
        import matplotlib.pyplot as plt
        
        ell = self.NlTT[:,0]
        f, axs = plt.subplots(1,2, sharex = True, sharey = True)
        color = ['c','b','k','y','g','r','gold','sienna','coral','navy']
        if list_freqs is None:
            axs[0].loglog(ell,self.NlTT[:,1]*2*np.pi/ell/(ell+1), label = "Primary channel")
            axs[1].loglog(ell,self.NlEE[:,1]*2*np.pi/ell/(ell+1), label = "Primary channel")
        else :
            k = 0
            for fr in list_freqs:
                i = self.freqs.index(fr)
                axs[0].loglog(ell,self.NlTT[:,2+i]*2*np.pi/ell/(ell+1), c = color[k], label = "{:d} GHz".format(fr))
                axs[1].loglog(ell,self.NlEE[:,2+i]*2*np.pi/ell/(ell+1), c = color[k], label = "{:d} GHz".format(fr))
                k+=1
        axs[0].set_title("N($\ell$) Temperature")
        axs[1].set_title("N($\ell$) Polarization")
        axs[0].set_ylabel("N($\ell$) [$\mu K^2$-SR] " )
        axs[1].set_ylabel("N($\ell$) [$\mu K^2$-SR] " )
        axs[0].set_xlabel("$\ell$")
        axs[1].set_xlabel("$\ell$")
        axs[0].set_ylim(1e-6,10)
        axs[0].legend(loc = 'upper right')
        axs[1].legend(loc = 'upper right')
        axs[0].set_xlim(10, int(ell[-1]))
        plt.suptitle("Noise power spectra for {}".format(self.name))
        plt.show()


class Setup():
    """ Class that summarizes the setup used for the Fisher Forecast
        - param_list : lsit of paramters used for Fisher analysis
        - noise_root : directory where noise data are written
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
    def __init__(self,param_list,noise_root,data_root, use_BBN, mass_neutrinos = 60., l_boost = 1, boost = 1):
        """ Initialization """
        self.param_list = param_list
        self.noise_root = noise_root
        self.data_root = data_root
        self.use_BBN = use_BBN
        self.mass_neutrinos = mass_neutrinos
        self.l_boost = l_boost
        self.boost = boost
        
    def get_list_experiments(self, list_experiments):
        """ Method to get the list of combined experiment in the setup. Separated from the rest of the __init__ method to make parsing easier"""
        self.list_experiments = list_experiments
        self.lmax = max([max(expe.lmax_T,expe.lmax_P) for expe in list_experiments])
    
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
    print(param_list)
    try:
        assert set(param_list).issubset(['H0','A_s','109A_s','ln10A_s','N_eff','Y_He','ombh2','omch2','theta_MC','n_s','tau'])
    except AssertionError:
        print("Some parameters are invalid. Check you ini file")
    noise_root = config.get(setup_sec,'noise_root')
    data_root = config.get(setup_sec,'data_root')
    use_BBN = config.getfloat(setup_sec,'use_BBN')
    if config.has_option(setup_sec,'mass_neutrinos'):
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
    
    setup = Setup(param_list,noise_root,data_root, use_BBN, mass_neutrinos, l_boost, boost)
            
    list_experiments_included = []
       
    ## READING EXPERIMENTS CONFIGURATIONS ##
    
    for experiment in sections[1:]:
        atmo_noise = False
        SO = False
        name = config.get(experiment,'name')
        fsky = config.getfloat(experiment,'fsky')
        include = config.getboolean(experiment,'include')
        include_P = config.getboolean(experiment,'polarization')
        include_lens = config.getboolean(experiment,'lensing')
        include_rayleigh = config.getboolean(experiment,'rayleigh')
        lmin = config.getint(experiment,'lmin')
        lmax_T = config.getint(experiment,'lmax_T')
        lmax_P = config.getint(experiment,'lmax_P')
        freqs = [ int(fr) for fr in config.get(experiment,'freqs').split(',')]      
        if config.has_option(experiment,'sensitivity'):
            sensitivity = config.getint(experiment,'sensitivity')
            SO = True
        else:
            if config.has_option(experiment,'noise_pixel_T'):
                noise_pix_T = [float(noise) for noise in config.get(experiment,'noise_pixel_T').split(',')]
            if config.has_option(experiment,'noise_pixel_P'):
                noise_pix_P = [float(noise) for noise in config.get(experiment,'noise_pixel_P').split(',')]
            else:
                noise_pix_P = None
            if config.has_option(experiment,'beam_FWHM'):
                beam_FWHM = [float(beam) for beam in config.get(experiment,'beam_FWHM').split(',')]
        #Parameters for atmospheric noise, all must be specified !
        if config.has_option(experiment,'alpha_temp'):
            atmo_noise = True
            alpha_temp = [float(alpha) for alpha in config.get(experiment,'alpha_temp').split(',')]
            if config.has_option(experiment,'alpha_pol'):
                alpha_pol = [float(alpha) for alpha in config.get(experiment,'alpha_pol').split(',')]
            else : 
                print("alpha_pol must be defined when using atmospheric noise")
            if config.has_option(experiment,'ell_pivot_temp'):
                ell_pivot_temp = [int(ell) for ell in config.get(experiment,'ell_pivot_temp').split(',')]
            else : 
                print("ell_pivot_temp must be defined when using atmospheric noise")
            if config.has_option(experiment,'ell_pivot_pol'):
                ell_pivot_pol = [int(ell) for ell in config.get(experiment,'ell_pivot_pol').split(',')]         
            else : 
                print("ell_pivot_pol_pol must be defined when using atmospheric noise")
            if config.has_option(experiment,'c_atmo_temp'):
                c_atmo_temp = [float(c_atmo) for c_atmo in config.get(experiment,'c_atmo_temp').split(',')]
            else : 
                print("c_atmo_temp must be defined when using atmospheric noise"  ) 
                    
        expe = Experiment(name,include,include_P,include_lens,include_rayleigh,freqs,fsky,lmax_T,lmax_P,lmin)
        if config.has_option(experiment,'lensing_estimators'):
            expe.lensing_estimators = config.getint(experiment,'lensing_estimators') 
        if config.has_option(experiment,'update'):
            expe.update = config.getboolean(experiment,'update')
        else:
            expe.update = True  
        expe.atmospheric = atmo_noise   
        if atmo_noise:
            expe.alpha_temp = alpha_temp
            expe.alpha_pol = alpha_pol
            expe.ell_pivot_temp = ell_pivot_temp
            expe.ell_pivot_pivot = ell_pivot_pol
            expe.c_atmo_temp = c_atmo_temp
            
        if SO: 
            expe.sensitivity = sensitivity
        else:
            expe.noise_pix_T = noise_pix_T
            expe.noise_pix_P = noise_pix_P
            expe.beam_FWHM = beam_FWHM
        if include:
            list_experiments_included.append(expe)
        
    
    return setup,list_experiments_included
    
    
if __name__ == "__main__":
    exp = Experiment('pla',True,True,True,[0,1,2],'../data/noise_spectra/template.dat',1.,100,100,2)
    exp.read_noise()

