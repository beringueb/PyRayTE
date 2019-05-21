#!/usr/bin/python3
"""Main module that drives the analysis"""

import numpy as np
import matplotlib.pyplot as plt
import experiment as ex
import CAMB_wrap as cw
import fisher as fi
import os

lmax_list = np.logspace(np.log10(31),4,30).astype('int')
#lmax_list = [50,100,250,500,1000,1500,2000]
update = False
fisher_matrices_root = '../data/fisher_matrices/PICO_lmax/'
for lmax in lmax_list[::-1]:
    for fsky in [0.75]:
        for rayleigh in ['True', 'False']:
            inifile = '../ini/PICO_lmax/PICO_{:d}_{:3.2f}_{:s}.ini'.format(lmax,fsky,rayleigh)
            setup = ex.ini_driver(inifile)
            setup.get_fiducial(fid_file = '../data/templates/fiducial_parameters.txt')
            cw.init_file(setup)
            fisher_list = []
            for experiment in setup.list_experiments : # loop on every combined experiments
                if update:
                    cw.parameter_files(setup, experiment)      
                    cw.compile_CAMB(experiment)
                    cw.run_CAMB(experiment)
                    update = False
                #experiment.plot()
                #print(experiment.NlPP)
                name = experiment.name
                fisher = fi.FisherMatrix(setup.param_list,experiment)
                fisher.get_fisher(setup)
                fisher_list.append(fisher)
            fisher_tot = fisher_list[0] + fisher_list[1] + fisher_list[2]
            fisher_tot.write_to_txt(fisher_matrices_root, name = 'PICO_lensed_pol_v1_{:d}'.format(lmax))
            update = False
            print([fish.fsky for fish in fisher_list])



                    
           
    
    
    




        

    
    
    
            



