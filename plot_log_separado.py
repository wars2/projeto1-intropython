# -*- coding: utf-8 -*-
"""
Created on Tue May 21 13:53:38 2019

@author: H_Mion
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

ds = pd.read_csv('DoubleMuRun2011A.csv')
invariant_mass_1 = ds['M']

no_bins = 300
# Let's calculate the logarithms of the masses and weighs.
inv_mass_log = np.log10(invariant_mass_1)
#print(inv_mass_log.min(),inv_mass_log.max())
#print("")
weights = []
for a in invariant_mass_1:
    weights.append(no_bins/np.log(10)/a)

# Let's plot the weighted histogram.
plt.title('The histogram of the invariant masses of two muons \n')
plt.hist(inv_mass_log, bins=no_bins, range=(-0.5,2.5))
plt.yscale('log')
plt.show()

plt.title('Picos Rô, ômega \n')
plt.hist(inv_mass_log, bins=no_bins, range=(-0.15,-0.05), weights=weights, color="red")
plt.yscale('log')
plt.show()

plt.title('Pico Phi \n')
plt.hist(inv_mass_log, bins=no_bins, range=(-0.09,0.09), weights=weights, color="red")
plt.yscale('log')
plt.show()

plt.title('Pico J/Psi \n')
plt.hist(inv_mass_log, bins=no_bins, range=(0.43,0.55), weights=weights, color="red")
plt.yscale('log')
plt.show()

plt.title('Pico Psi linha \n')
plt.hist(inv_mass_log, bins=no_bins, range=(0.51,0.62), weights=weights, color="red")
plt.yscale('log')
plt.show()

plt.title('Pico Upsilon \n')
plt.hist(inv_mass_log, bins=no_bins, range=(0.96,1.05), weights=weights, color="red")
plt.yscale('log')
plt.show()

plt.title('Pico Z \n')
plt.hist(inv_mass_log, bins=no_bins, range=(1.85,2.05), weights=weights, color="red")
plt.yscale('log')

# Naming the labels and the title.
plt.xlabel('log10(invariant mass) [log10(GeV)]')
plt.ylabel('Number of the events')

plt.show()