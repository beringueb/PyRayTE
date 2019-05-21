#!/usr/bin/python3
"""Main module that drives the analysis"""

import numpy as np
import matplotlib.pyplot as plt
import experiment as ex
import CAMB_wrap as cw
import fisher as fi
import os

update = True
fisher_matrices_root = '../data/fisher_matrices/'

for N_det in [1e4,1e5,1e6]:
    for fsky in [0.25,0.5,0.75]:
        for beam in [1,2,3,4]:
            inifile = '../ini/grid/grid_{:d}_{:3.2f}_{:d}.ini'.format(int(N_det),fsky,beam)
            print(inifile)
            setup = ex.ini_driver(inifile)
            setup.get_fiducial(fid_file = '../data/templates/fiducial_parameters.txt')
            cw.init_file(setup)
            fisher_list = []
            for experiment in setup.list_experiments : # loop on every combined experiments
                if update:
                    cw.parameter_files(setup, experiment)      
                    cw.compile_CAMB(experiment)
                    cw.run_CAMB(experiment)
                fisher = fi.FisherMatrix(setup.param_list,experiment)
                fisher.get_fisher(setup)
                fisher_list.append(fisher)
            fisher_tot = fisher_list[0]
            fisher_tot.write_to_txt(fisher_matrices_root, name = 'grid_{:d}_{:3.2f}_{:d}'.format(int(N_det),fsky,beam))
            update = False



                    
           
    
    
    




        

    
    
    
            



