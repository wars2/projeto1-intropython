# -*- coding: utf-8 -*-
"""
Created on Mon May 20 15:24:01 2019

@author: Hugo Borges
"""

'''
#parte 1
#codigo estraido de: https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html

import matplotlib.pyplot as plt
from scipy.misc import electrocardiogram
from scipy.signal import find_peaks

x = electrocardiogram()[2000:4000]
peaks, _ = find_peaks(x, height=0)
plt.plot(x)
plt.plot(peaks, x[peaks], "x")
plt.plot(np.zeros_like(x), "--", color="gray")
plt.show()
border = np.sin(np.linspace(0, 3 * np.pi, x.size))
peaks, _ = find_peaks(x, height=(-border, border))
plt.plot(x)
plt.plot(-border, "--", color="gray")
plt.plot(border, ":", color="gray")
plt.plot(peaks, x[peaks], "x")
plt.show()
peaks, _ = find_peaks(x, distance=150)
#np.diff(peaks)
plt.plot(x)
plt.plot(peaks, x[peaks], "x")
plt.show()
peaks, properties = find_peaks(x, prominence=(None, 0.6))
properties["prominences"].max()
plt.plot(x)
plt.plot(peaks, x[peaks], "x")
plt.show()
x = electrocardiogram()[17000:18000]
peaks, properties = find_peaks(x, prominence=1, width=20)
properties["prominences"], properties["widths"]
plt.plot(x)
plt.plot(peaks, x[peaks], "x")
plt.vlines(x=peaks, ymin=x[peaks] - properties["prominences"], ymax = x[peaks], color = "C1")
plt.hlines(y=properties["width_heights"], xmin=properties["left_ips"], xmax=properties["right_ips"], color = "C1")
plt.show()
'''

'''
# parte 2
# parte extraida da parte "Making the histogram" do notebbok

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ds = pd.read_csv('DoubleMuRun2011A.csv')

invariant_mass = np.sqrt(2*ds.pt1*ds.pt2*((np.cosh(ds.eta1-ds.eta2))-(np.cos(ds.phi1-ds.phi2))))

plt.hist(invariant_mass, bins=600, range=(0,110))
print("Valor mínimo = " + str(invariant_mass.min()) + "\nValor máximo = " + str(invariant_mass.max()))

# Let's name the axes and the title. Don't change these.
plt.xlabel('Invariant mass [GeV]')
plt.ylabel('Number of events')
plt.title('Histogram of invariant mass values of two muons. \n')
plt.xscale("log")
plt.yscale("log")
plt.show()
'''

#parte 3
# parte extraida da parte "The histogram of the whole data" do notebbok

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

ds = pd.read_csv('DoubleMuRun2011A.csv')
invariant_mass_1 = ds['M']

no_bins = 600
# Let's calculate the logarithms of the masses and weighs.
inv_mass_log = np.log10(invariant_mass_1)
weights = []
for a in invariant_mass_1:
    weights.append(no_bins/np.log(10)/a)

print("Valor mínimo = " + str(inv_mass_log.min()) + "\nValor máximo = " + str(inv_mass_log.max()))

# Let's plot the weighted histogram.
n_hist, bins_hist, patches = plt.hist(inv_mass_log, bins=no_bins, range=(-0.5,2.5), weights=weights, color="red")
plt.yscale('log')

# Naming the labels and the title.
plt.xlabel('log10(invariant mass) [log10(GeV)]')
plt.ylabel('Number of the events')
plt.title('The histogram of the invariant masses of two muons \n')
plt.show()


peaks, _ = find_peaks(n_hist, distance=150)
np.diff(peaks)
plt.plot(x)
plt.plot(peaks, x[peaks], "x")
plt.show()


