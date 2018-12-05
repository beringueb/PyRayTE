#!/usr/bin/python3
"""Main module that drives the analysis"""

import numpy as np
import matplotlib.pyplot as plt
import experiment as ex
import CAMB_wrap as cw
import fisher as fi
import pickle
import os
update = True
inifile = ''
fisher_matrices_root = '../data/fisher_matrices/'

setup = ex.inidriver(inifile)
setup.get_fiducial(fid_file = '../data/templates/fiducial_parameters.txt')
cw.init_file(setup)
fisher_list = []
for experiment in setup.list_experiments : # loop on every combined experiments
    if experiment.update:
        cw.parameter_files(setup, experiment)      
        cw.compile_CAMB(experiment)
        cw.run_CAMB(experiment)
    name = experiment.name
    fisher = fi.FisherMatrix(setup.param_list,experiment)
    fisher.get_fisher(setup)
    fisher.write_to_txt(os.path.join(fisher_matrices_root,'{}.fish'.format(name)))
    fisher_list.append(fisher)

if __name__ == "__main__": 

                    
           
    
    
    




        

    
    
    
            



