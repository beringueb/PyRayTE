#!/bin/usr/python3
""" Module CAMB_wrap that provides a wrapper for using CAMB. It will generate one .ini file for the fiducial model and two additional .ini files, one for each direction of a step. 
    It requires the python module camb only for setting background cosmology. The reason camb is not called via pycamb is that pycamb does not include rayleigh scattering.
"""

import numpy as np
import os 
from shutil import copyfile
import time

try:
    import subprocess 
except ImportError:
    print("Module suprocess needs to be installed in order to call CAMB")
try:
    import camb
except ImportError:
    print("pycamb (module camb) needs to be installed in order to set background cosmology")
    
CAMB_ROOT = "/Users/benjamin_brng/Documents/Cambridge/PhD/Rayleigh/Code/CAMB/"  # Location of the CAMB folder.

def init_file(setup):
    """ Function to generate the .ini file for the fiducial cosmology.
        - setup : setup object defined in 
    """ ### TO COMPLETE
    parameter_list = setup.param_list
    pars_fid = camb.CAMBparams()
    
    if 'H0' in parameter_list:
        H0_tmp = setup.fiducial['H0']
        cosmomc_theta_tmp = None
    else:
        H0_tmp = None
        cosmomc_theta_tmp=setup.fiducial['theta_MC']
    if 'Y_He' in parameter_list:
        YHe_tmp = setup.fiducial['Y_He']
    elif setup.use_BBN == 1:
        YHe_tmp = None #Y_p is then fixed by BBN consistency
    else:
        YHe_tmp = setup.use_BBN #Y_p is fixed but not by BBN consistency
     
    if 'ln10A_s' in parameter_list:
        As = np.exp(setup.fiducial['ln10A_s'])*1.E-10
    elif '109A_s' in parameter_list:
        As = setup.fiducial['109A_s']*1.E-9
    else:
        print('Due to numerical instabilities, it is not advised to directly use A_s as a parameter, use ln10A_s or 109A_s')
        As = setup.fiducial['A_s']
    
    #Setting cosmology in CAMB    
    mass_nu = setup.fiducial['mass_neutrinos'] if 'mass_neutrinos' in setup.param_list else setup.mass_neutrinos
    pars_fid.set_cosmology(H0=H0_tmp,cosmomc_theta=cosmomc_theta_tmp,ombh2=setup.fiducial['ombh2'],omch2=setup.fiducial['omch2'],mnu = mass_nu,nnu=setup.fiducial['N_eff'],tau=setup.fiducial['tau'],YHe=YHe_tmp)
    #Initializing parameter file for CAMB
    param_file = os.path.join(CAMB_ROOT,"params.ini")
    new_file = ''
    with open(param_file,'r') as f :
        for line in f :
            if 'hubble' in line and not('#' in line):
                line1,line2 = line.split('=')
                line2 = '=  {} \n'.format(pars_fid.H0)
                line = line1 + line2
            if 'ombh2' in line and not('#' in line):
                line1,line2 = line.split('=')
                line2 = '=  {} \n'.format(setup.fiducial['ombh2'])
                line = line1 + line2
            if 'omch2' in line and not('#' in line):
                line1,line2 = line.split('=')
                line2 = '=  {} \n'.format(setup.fiducial['omch2'])
                line = line1 + line2
            if 'omnuh2' in line and not('#' in line):
                line1,line2 = line.split('=')
                line2 = '=  {} \n'.format(pars_fid.omnuh2 )
                line = line1 + line2
            if 'omk' in line and not('#' in line):
                line1,line2 = line.split('=')
                line2 = '=  {} \n'.format(pars_fid.omk)
                line = line1 + line2
            if 'massless_neutrinos' in line and not('#' in line):
                line1,line2 = line.split('=')
                line2 = '=  {} \n'.format(setup.fiducial['N_eff'] - pars_fid.num_nu_massive)
                line = line1 + line2
            if 'massive_neutrinos' in line and not('#' in line):
                line1,line2 = line.split('=')
                line2 = '=  {} \n'.format(pars_fid.num_nu_massive)
                line = line1 + line2
            if 'l_max_scalar' in line and not('#' in line):
                line1,line2 = line.split('=')
                line2 = '=  {} \n'.format(setup.lmax+500) #Need to check up to which l to go 
                line = line1 + line2
            if 'scalar_amp' in line and not('#' in line):
                line1,line2 = line.split('=')
                line2 = '=  {} \n'.format(As)
                line = line1 + line2
            if 'scalar_spectral_index' in line and not('#' in line):
                line1,line2 = line.split('=')
                line2 = '=  {} \n'.format(setup.fiducial['n_s']) 
                line = line1 + line2
            if 're_optical_depth' in line and not('#' in line):
                line1,line2 = line.split('=')
                line2 = '=  {} \n'.format(pars_fid.Reion.optical_depth)
                line = line1 + line2
            if 'helium_fraction' in line and not('#' in line):
                line1,line2 = line.split('=')
                line2 = '=  {} \n'.format(pars_fid.YHe)
                line = line1 + line2
            if 'output_root' in line and not('#' in line):
                line1,line2 = line.split('=')
                line2 = '=  {} \n'.format(setup.data_root)
                line = line1 + line2
            if 'l_accuracy_boost' in line and not('#' in line):
                line1,line2 = line.split('=')
                line2 = '=  {}  \n'.format(setup.l_boost)
                line = line1 + line2
            if 'accuracy_boost' in line and not('#' in line) and not('l_accuracy_boost' in line):
                line1,line2 = line.split('=')
                line2 = '=  {} \n'.format(setup.boost)
                line = line1 + line2
            if 'lensed_output_file' in line and not('#' in line):
                line1,line2 = line.split('=')
                line2 = '=  {} \n'.format("scalCls_lensed") 
                line = line1 + line2
            if 'lens_potential_output_file' in line and not('#' in line):
                line1,line2 = line.split('=')
                line2 = '=  {} \n'.format("scalCls_lenspot.dat")
                line = line1 + line2
            new_file += line
    param_file_init = os.path.join(CAMB_ROOT,"params_init.ini")
    with open(param_file_init, "w") as f:
        f.write(new_file)
    print('Initial parameter file written at {}'.format(CAMB_ROOT))

def parameter_files(setup, experiment):
    """ Function that creates all the .ini files to compute the derivatives and the fiducial power spectra. These files are located in a directory named 'Experiment_inifiles'. This directory is emptied and populated again. Then CAMB is run for each file in this directory.
        - setup : 
        - experiment : 
    """ # TO COMPLETE
    # First check if the directory exists, if not create it, if yes, empty it
    if not os.path.exists(os.path.join(CAMB_ROOT,"{}_inifiles".format(experiment.name))):
        os.mkdir(os.path.join(CAMB_ROOT,"{}_inifiles".format(experiment.name))) 
    else : 
        list_file = os.listdir(os.path.join(CAMB_ROOT,"{}_inifiles".format(experiment.name)))
        for f in list_file:
            os.unlink(os.path.join(CAMB_ROOT,"{}_inifiles".format(experiment.name),f))
        try:
            assert not os.listdir(os.path.join(CAMB_ROOT,"{}_inifiles".format(experiment.name)))
        except AssertionError:
            print("ini files directory for {} is not empty, something is wrong !".format(experiment.name))
     
    # Do the same for the data files
    if not os.path.exists(os.path.join(setup.data_root,experiment.name)):
        os.mkdir(os.path.join(setup.data_root,experiment.name)) 
    else : 
        list_file = os.listdir(os.path.join(setup.data_root,experiment.name))
        for f in list_file:
            os.unlink(os.path.join(setup.data_root,experiment.name,f))
        try:
            assert not os.listdir(os.path.join(setup.data_root,experiment.name))
        except AssertionError:
            print("Data directory for {} is not empty, something is wrong !".format(experiment.name))
    
    # Get a copy of the fiducial .ini file and change the output root to get a copy of the fiducial power spectra in the directory
    copyfile(os.path.join(CAMB_ROOT,"params_init.ini"),os.path.join(CAMB_ROOT,"{}_inifiles".format(experiment.name),"params_fiducial.ini"))
    newText = ''
    with open(os.path.join(CAMB_ROOT,"{}_inifiles".format(experiment.name),"params_fiducial.ini"),'r') as f :
        for line in f : 
            if 'output_root' in line and not('#' in line):
                line1,line2 = line.split('=')
                line2 = '=  {} \n'.format(os.path.join(setup.data_root,experiment.name,"fiducial"))
                line = line1 + line2
            newText += line
    with open(os.path.join(CAMB_ROOT,"{}_inifiles".format(experiment.name),"params_fiducial.ini"),'w') as f :
        f.write(newText)

   
    # Then for each parameter in parameter_list, generate two files, one in each direction to get the derivatives
    parameter_list = setup.param_list
    for param in parameter_list:
        for direction in ['left','right']:
            pars_temp = camb.CAMBparams()
            temp = setup.fiducial.copy() # take a copy of the fiducial value of the parameter
            if direction is 'left':
                temp[param] = setup.fiducial[param] - setup.step[param]
            else : 
                temp[param] = setup.fiducial[param] + setup.step[param]
            if 'H0' in parameter_list:
                H0_tmp = temp['H0']
                cosmomc_theta_tmp = None
            else:
                H0_tmp = None
                cosmomc_theta_tmp=temp['theta_MC']
            if 'Y_He' in parameter_list:
                YHe_tmp = temp['Y_He']
            elif setup.use_BBN == 1:
                YHe_tmp = None #Y_p is then fixed by BBN consistency
            else:
                YHe_tmp = setup.use_BBN #Y_p is fixed but not by BBN consistency
     
            if 'ln10A_s' in parameter_list:
                As = np.exp(temp['ln10A_s'])*1E-10
            elif '109A_s' in parameter_list:
                As = temp['109A_s']*1E-9
            else:
                As = temp['A_s']
            
            mass_nu = temp['mass_neutrinos'] if 'mass_neutrinos' in setup.param_list else setup.mass_neutrinos
            
            pars_temp.set_cosmology(H0=H0_tmp,cosmomc_theta=cosmomc_theta_tmp,ombh2=temp['ombh2'],omch2=temp['omch2'],mnu = mass_nu,nnu=temp['N_eff'],tau=temp['tau'],YHe=YHe_tmp)
            param_file = os.path.join(CAMB_ROOT,"{}_inifiles".format(experiment.name),"params_fiducial.ini")
            new_file = ''
            with open(param_file,'r') as f :
                for line in f :
                    if 'hubble' in line and not('#' in line):
                            line1,line2 = line.split('=')
                            line2 = '=  {} \n'.format(pars_temp.H0)
                            line = line1 + line2
                    if 'ombh2' in line and not('#' in line):
                        line1,line2 = line.split('=')
                        line2 = '=  {} \n'.format(temp['ombh2']) # .format(pars_temp.ombh2)
                        line = line1 + line2
                    if 'omch2' in line and not('#' in line):
                        line1,line2 = line.split('=')
                        line2 = '=  {} \n'.format(temp['omch2']) #.format(pars_temp.omch2) #
                        line = line1 + line2#
                    if 'omnuh2' in line and not('#' in line):
                        line1,line2 = line.split('=')
                        line2 = '=  {} \n'.format(pars_temp.omnuh2)
                        #line = line1 + line2
                    if 'massless_neutrinos' in line and not('#' in line):
                        line1,line2 = line.split('=')
                        line2 = '=  {} \n'.format(temp['N_eff'] - pars_temp.num_nu_massive)
                        line = line1 + line2
                    if 'massive_neutrinos' in line and not('#' in line):
                        line1,line2 = line.split('=')
                        line2 = '=  {} \n'.format(pars_temp.num_nu_massive)
                        line = line1 + line2
                    if 'scalar_amp' in line and not('#' in line):
                        line1,line2 = line.split('=')
                        line2 = '=  {} \n'.format(As)
                        line = line1 + line2
                    if 'scalar_spectral_index' in line and not('#' in line):
                        line1,line2 = line.split('=')
                        line2 = '=  {} \n'.format(temp['n_s']) 
                        line = line1 + line2
                    if 're_optical_depth' in line and not('#' in line):
                        line1,line2 = line.split('=')
                        line2 = '=  {} \n'.format(pars_temp.Reion.optical_depth)
                        line = line1 + line2
                    if 'helium_fraction' in line and not('#' in line):
                        line1,line2 = line.split('=')
                        line2 = '=  {} \n'.format(pars_temp.YHe)
                        line = line1 + line2
                    if 'output_root' in line and not('#' in line):
                        line1,line2 = line.split('=')
                        line2 = '=  {} \n'.format(os.path.join(setup.data_root,experiment.name,"{}_{}".format(param,direction)))
                        line = line1 + line2

                    new_file += line
            param_file_out = os.path.join(CAMB_ROOT,"{}_inifiles".format(experiment.name),"params_{}_{}.ini".format(param,direction))
            with open(param_file_out, "w") as f:
                f.write(new_file)
    print('Initial parameter files written for {}'.format(experiment.name))
    

def compile_CAMB(experiment):
    """ Function to compile CAMB for the correct frequencies
        - experiment: 
    """ # TO COMPLETE
    
    # Information about the frequencies used by CAMB are located in modules.f90
    file_module = os.path.join(CAMB_ROOT,"modules.f90")
    newText = ''
    freqs_write = experiment.freqs.copy()
    freqs_write.insert(0,0) # insert freq 0 GHz to get the no Rayleigh case
    with open(file_module,'r') as f:
        for line in f:
            if 'num_cmb_freq =' in line:
                line1,line2 = line.split('=')
                line2 = '=  {:d} \n'.format(len(experiment.freqs)+1)
                line = line1 + line2
            if 'phot_freqs(1:' in line :
                line1,line2 = line.split('(1:')
                line2 = '(1:{:d}) = {} \n'.format(len(experiment.freqs)+1,str(freqs_write)) 
                line = line1 + line2
            newText += line
    with open(file_module,'w') as f:
        f.write(newText)
    print("modules.f90 file updated, compiling")
    subprocess.call("make clean", shell = True, cwd = CAMB_ROOT, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
    subprocess.call("make", shell = True, cwd = CAMB_ROOT, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
    print("Done ! ")
    

def run_CAMB(experiment):
    """ Function to run CAMB on every file located at 'experiment_inifiles' 
        - experiment : 
    """ # TO COMPLETE
    
    param_files_list = os.listdir(os.path.join(CAMB_ROOT,"{}_inifiles".format(experiment.name))) 
    n_file = len(param_files_list)
    print("Running CAMB on {:d} files ... ".format(n_file))
    time_start = time.time()
    i = 1.
    for f in param_files_list : # loop over all the .ini files
        file_name = os.path.join(CAMB_ROOT,"{}_inifiles".format(experiment.name),f)
        subprocess.call("./camb {}".format(file_name), shell = True, cwd = CAMB_ROOT, stdout = subprocess.DEVNULL) # CHANGE HYREC to output in stdout and not stderr
        time_tmp = time.time()
        ETA = (n_file - i)*(time_tmp - time_start) / i
        print("{:3.1f}% done, ETA : {:2.0f} min {:2.0f} secs".format(i/n_file * 100, ETA // 60, ETA % 60), end = "\r" )
        i += 1
    print("Done in {:2.0f} min {:2.0f} secs".format((time_tmp - time_start) // 60, (time_tmp - time_start) % 60))
    

    
        
        
         
             
            
            
                
            
            
            
            
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
    



