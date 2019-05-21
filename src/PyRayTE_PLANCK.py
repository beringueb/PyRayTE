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

inifile = '../ini/PLANCK.ini'
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
    fisher.write_to_txt(fisher_matrices_root, name = '{}'.format(name))
    fisher_list.append(fisher)
    
    
fisher_tot = fisher_list[0] + fisher_list[1] 
fisher_tot.write_to_txt(fisher_matrices_root, name = 'PLANCK_only')

print([fish.fsky for fish in fisher_list])
print(fi.print_errors([fisher_tot]))


                    
           
    
    
    




        

    
    
    
            



