#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 18:47:26 2020

@author: mikkel
"""
import os
if "StatMechGlass" in os.getcwd():
    import stat_mech_glass as smg
else:
    # This is how you want to import the package
    from StatMechGlass import stat_mech_glass as smg

# Na.csv and Na_Tg.csv are located in Data/SiO2. This command will build
# and save the Na-Si interaction parameters
smg.smg_binary_par("Si", "Na", it=10)

# Na.csv is located in Data/SiB. This command will build and save the Si/B 
# interaction parameter, using Na2O-SiO2-B2O3 data (any modifier would do)
smg.smg_ternary_par(["Si","B"], "Na", it=1) 

# The glass composition is defined using a python dictionary. 
# This would correspond to a 25Na2O-25B2O3-25SiO2 glass
glass_comp = {"Si": 25, "B": 25, "Na": 25}
# Results are stored in a variable for any desired use 
results = smg.smg_structure(glass_comp, 700)
print(results)

# This will make a plot of the structures in the 25B2O3-25SiO2 glass system
# as a function of Na2O content 
smg.smg_plot(glass_comp, "Na", 700, plt_save = False)