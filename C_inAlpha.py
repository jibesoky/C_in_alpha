# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import os
from numpy.polynomial import polynomial as P

files = os.listdir()
my_files = []
my_data = {}

# The C_inAlpa_%.txt files were obtained from Thermocalc paraequilibrium calculations between fcc and bcc phases.

for i in files:
    if i.startswith('C_inAlpha_'):
        filename = str(i)
        my_files.append(filename)
        my_data[filename] = np.loadtxt(filename)

my_data_sort = {}
poly_adj = {}
C_in_alpha_fit ={}
temp=np.linspace(73,873,100)
for i in my_files:
    my_data_sort[i] = my_data[i][np.argsort(my_data[i][:,0])]
    poly_adj[i] = np.polyfit(my_data_sort[i][:,0]-273.15, my_data_sort[i][:,4],6)
    C_in_alpha_fit[i] =  np.polyval(poly_adj[i],temp)
    plt.plot(my_data_sort[i][:,0]-273.15, my_data_sort[i][:,4],'o',label=i)
    plt.plot(temp, np.polyval(poly_adj[i],temp),label='Fit of ' + i)

def C_in_alpha(T): #Mole fraction of C in alpha using thermocalc under para equilibrium
    C= 1.4734491E-20*T**6 + 3.9638142E-17*T**5 - 1.1293268E-13*T**4 + 6.8406210E-11*T**3 - 9.3489472E-09*T**2 + 6.1810195E-07*T - 6.3920771E-06
    return abs(C)

MF_to_WP=np.array([4.88194823e-05,2.14076779e+01,1.69714954e+01,1.16633462e+01,1.70719962e+01])
def mf_to_wp(MF):
    return P.polyval(MF,MF_to_WP)  

C=C_in_alpha(temp)
# Cwp=mf_to_wp(C)

plt.plot(temp,C,'y', label='fit from master fitter')
plt.legend()


