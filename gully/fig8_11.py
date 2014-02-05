"""
Cosmology Regression Example
----------------------------
Figure 8.11

A Gaussian process regression analysis of the simulated supernova sample used
in figure 8.2. This uses a squared-exponential covariance model, with bandwidth
learned through cross-validation.
"""
# Author: Jake VanderPlas
# License: BSD
#   The figure produced by this code is published in the textbook
#   "Statistics, Data Mining, and Machine Learning in Astronomy" (2013)
#   For more information, see http://astroML.github.com
#   To report a bug or issue, use the following forum:
#    https://groups.google.com/forum/#!forum/astroml-general

#Edited on Feb 4, 2014 by gully for SSC383D
#For exercise 2 part (A)
#Pseudocode:
#1) Assemble a training dataset and a cross validation dataset
#2) Fit them with a variety of h smoothing kernels
#3) Calculate the mean squared prediction error for each one


import numpy as np
from matplotlib import pyplot as plt
from scipy import ndimage as ndimage

#----------------------------------------------------------------------
# This function adjusts matplotlib settings for a uniform feel in the textbook.
# Note that with usetex=True, fonts are rendered with LaTeX.  This may
# result in an error if LaTeX is not installed on your system.  In that case,
# you can set usetex to False.
from astroML.plotting import setup_text_plots
setup_text_plots(fontsize=16, usetex=False)

#------------------------------------------------------------
# Generate data
x=np.linspace(0.0,1.0, 500)
amp=1.0
p1=0.1
p2=1.0
nz1=0.05
nz2=0.35

#                   high noise      low noise
#Wiggly function    11              12
#smooth function    21              22

y11= np.sin(2.0*np.pi*x/p1) + np.random.normal(loc=0.0, scale=nz2, size=500)
y12= np.sin(2.0*np.pi*x/p1) + np.random.normal(loc=0.0, scale=nz1, size=500)
y21= np.sin(2.0*np.pi*x/p2) + np.random.normal(loc=0.0, scale=nz2, size=500)
y22= np.sin(2.0*np.pi*x/p2) + np.random.normal(loc=0.0, scale=nz1, size=500)

#cross validation:
w_arr=np.linspace(0, 99, num=100)+3
cv_arr11=w_arr*0.0
cv_arr22=w_arr*0.0
for i in range(0, 100, 1):
    sm_w=w_arr[i]
    inds1=range(0, 499, 2)
    inds2=range(1, 500, 2)
    y11_s = ndimage.filters.gaussian_filter1d(y11[inds1], sm_w)
    y11_cv= y11[inds2] - y11_s
    cv_mse= np.mean(y11_cv**2.0)
    cv_arr11[i]=cv_mse

for i in range(0, 100, 1):
    sm_w=w_arr[i]
    inds1=range(0, 499, 2)
    inds2=range(1, 500, 2)
    y22_s = ndimage.filters.gaussian_filter1d(y22[inds1], sm_w)
    y22_cv= y22[inds2] - y22_s
    cv_mse= np.mean(y22_cv**2.0)
    cv_arr22[i]=cv_mse
    
#------------------------------------------------------------
# Plot 
fig = plt.figure(figsize=(7, 7))
fig.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)

ax11 = fig.add_subplot(221)
ax11.plot(x[inds1], y11[inds1], '.', color='#DDDDDD')
ax11.plot(x[inds1], y11[inds2], '.', color='#0000FF')
ax11.plot(x[inds1], y11_s, '-k')

ax11cv = fig.add_subplot(222)
ax11cv.plot(w_arr, cv_arr11, '-r')

ax22 = fig.add_subplot(223)
ax22.plot(x[inds1], y22[inds1], '.', color='#DDDDDD')
ax22.plot(x[inds1], y22[inds2], '.', color='#0000FF')
ax22.plot(x[inds1], y22_s, '-k')

ax22cv = fig.add_subplot(224)
ax22cv.plot(w_arr, cv_arr22, '-r')
#ax22.plot(x, y22, '-r')
#ax22.plot(x, y22, '-r')
#ax12.plot(x, y12, '-r')

#ax21 = fig.add_subplot(223)
#ax21.plot(x, y21, '-r')

#ax22 = fig.add_subplot(224)
#ax22.plot(x, y22, '-r')

#ax.plot(z_fit, y_pred, '-k')
#ax.fill_between(z_fit, y_pred - 1.96 * sigma, y_pred + 1.96 * sigma,
#                alpha=0.2, color='b', label='95% confidence interval')

ax11.set_ylim(-1.5, 1.5)
ax12.set_ylim(-1.5, 1.5)
ax21.set_ylim(-1.5, 1.5)
ax22.set_ylim(-1.5, 1.5)

plt.show()
