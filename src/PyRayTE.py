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

for sensi in [0,1,2]:
    for fsky in [0.1,0.2,0.4]:
        inifile = '../ini/SO/SO_{:d}_{:3.1f}.ini'.format(sensi,fsky)
        setup = ex.ini_driver(inifile)
        setup.get_fiducial(fid_file = '../data/templates/fiducial_parameters.txt')
        cw.init_file(setup)
        fisher_list = []
        for experiment in setup.list_experiments : # loop on every combined experiments
            if update:
                cw.parameter_files(setup, experiment)      
                cw.compile_CAMB(experiment)
                cw.run_CAMB(experiment)
            #experiment.plot()
            #print(experiment.NlPP)
            name = experiment.name
            fisher = fi.FisherMatrix(setup.param_list,experiment)
            fisher.get_fisher(setup)
            fisher_list.append(fisher)
        fisher_tot = fisher_list[0] + fisher_list[1] + fisher_list[2]
        fisher_tot.write_to_txt(fisher_matrices_root, name = 'SO+PLANCK_{:d}_{:3.1f}'.format(sensi,fsky))
        update = False
        print([fish.fsky for fish in fisher_list])



                    
           
    
    
    




        

    
    
    
            



