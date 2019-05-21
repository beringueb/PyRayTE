#!/usr/bin/python3
"""Main module that drives the analysis"""

import numpy as np
import matplotlib.pyplot as plt
import experiment as ex
import CAMB_wrap as cw
import fisher as fi
import os

update = False
fisher_matrices_root = '../data/fisher_matrices/'

inifile = '../ini/SO_test.ini'
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
fisher_tot = fisher_list[0] 
fisher_tot.write_to_txt(fisher_matrices_root, name = 'SO_test')
update = False
print([fish.fsky for fish in fisher_list])



                    
           
    
    
    




        

    
    
    
            



