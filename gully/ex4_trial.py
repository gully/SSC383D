# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 10:38:40 2014

@author: gully
"""

import numpy as np
from matplotlib import pyplot as plt
import csv
import os as os
import statsmodels.sandbox.gam as gam


#----------------------------------------------------------------------
# This function adjusts matplotlib settings for a uniform feel in the textbook.
# Note that with usetex=True, fonts are rendered with LaTeX.  This may
# result in an error if LaTeX is not installed on your system.  In that case,
# you can set usetex to False.
from astroML.plotting import setup_text_plots
setup_text_plots(fontsize=12, usetex=False) 
#MGS Note: change to 'True' if you have Latex issues


#Read in the data
#-----------------------------------------
dir1='/Users/gully/astroML/SSC383D/SSC383D/Chapter04/examples/'
os.chdir(dir1)
fn='air.csv'

ifile  = open(fn, "rb")
reader = csv.reader(ifile)

dat=np.zeros((6, 112))

rownum = 0
for row in reader:
    # Save header row.
    if rownum == 0:
        header = row
    else:
        colnum = 0
        for col in row:
            dat[colnum,rownum]=col
            colnum += 1
            
    rownum += 1

ifile.close()
#-----------------------------------------


#Wrangle the data
#-----------------------------------------
ndat=dat[:,1:]
ozone=ndat[0,:]
solarR=ndat[1,:]
wind=ndat[2,:]
temp=ndat[3,:]
#-----------------------------------------


#Plot it
#-----------------------------------------
plt.figure()               
ax=plt.axes()               

print plt.get_fignums()
print plt.get_backend()

plt.axes(axisbg='#eeeeee')
#plt.plot(temp, ozone, '.',color='#0000FF')

#plt.plot(solarR, ozone, '.',color='#0000FF')

#plt.xlabel(r'Temperature (F)', fontsize=16)
#plt.xlabel(r'Solar Insolation (Ly)', fontsize=16)
#plt.xlabel(r'Wind (mph)', fontsize=16)
plt.xlabel('Fitted Ozone Conc. (ppb)', fontsize=16)
plt.ylabel('Ozone Conc. (ppb)', fontsize=16)
plt.title(r'Ozone regression example', fontsize=20)

#plt.savefig('ex4_fig.eps', format='eps', facecolor='#eeeeee',
#            edgecolor='#eeeeee', )
#plt.rcdefaults()
#-----------------------------------------

#Fit an additive model
#ycoeff=np.polyfit(temp, ozone, 1)
ycoeff1=np.polyfit(solarR, ozone, 1)
print ycoeff1
p1=np.poly1d(ycoeff1)

ycoeff2=np.polyfit(wind, ozone, 1)
print ycoeff2
p2=np.poly1d(ycoeff2)

ycoeff3=np.polyfit(temp, ozone, 1)
print ycoeff3
p3=np.poly1d(ycoeff3)

yfit=p1(solarR)+ p2(wind)+ p3(temp)-88.0
plt.plot(yfit, ozone, '.',color='#00FF00')
plt.plot([0, 250], [0, 250], '-',color='#FF0000')

plt.show()
print 'the end'