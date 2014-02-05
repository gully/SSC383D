"""
Cosmology Regression Example
----------------------------
Figure 8.2

Various regression fits to the distance modulus vs. redshift relation for a
simulated set of 100 supernovas, selected from a distribution
:math:`p(z) \propto (z/z_0)^2 \exp[(z/z_0)^{1.5}]` with :math:`z_0 = 0.3`.
Gaussian basis functions have 15 Gaussians evenly spaced between z = 0 and 2,
with widths of 0.14. Kernel regression uses a Gaussian kernel with width 0.1.
"""
# Author: Jake VanderPlas
# License: BSD
#   The figure produced by this code is published in the textbook
#   "Statistics, Data Mining, and Machine Learning in Astronomy" (2013)
#   For more information, see http://astroML.github.com
#   To report a bug or issue, use the following forum:
#    https://groups.google.com/forum/#!forum/astroml-general
import numpy as np
from matplotlib import pyplot as plt

from astroML.cosmology import Cosmology
from astroML.datasets import generate_mu_z
from astroML.linear_model import NadarayaWatson

#----------------------------------------------------------------------
# This function adjusts matplotlib settings for a uniform feel in the textbook.
# Note that with usetex=True, fonts are rendered with LaTeX.  This may
# result in an error if LaTeX is not installed on your system.  In that case,
# you can set usetex to False.
from astroML.plotting import setup_text_plots
setup_text_plots(fontsize=12, usetex=False)

#------------------------------------------------------------
# Generate data
z_sample, mu_sample, dmu = generate_mu_z(100, random_state=0)

cosmo = Cosmology()
z = np.linspace(0.01, 2, 1000)
mu_true = np.asarray(map(cosmo.mu, z))

#------------------------------------------------------------
# Define our classifiers
basis_mu = np.linspace(0, 2, 15)[:, None]
basis_sigma = 3 * (basis_mu[1] - basis_mu[0])

subplots = [121, 122]


# number of constraints of the model.  Because
# Nadaraya-watson is just a weighted mean, it has only one constraint

#------------------------------------------------------------
# Plot the results
NN=20
h_arr=np.linspace(0.03, 0.5, num=NN)
crossval=h_arr*0.0

fig = plt.figure(figsize=(9, 5))

fig.subplots_adjust(left=0.2, right=0.95,
                bottom=0.15, top=0.95,
                hspace=0.1, wspace=0.3)
ax = fig.add_subplot(121)
ax2= fig.add_subplot(122)

for i in range(0, NN,1):

    subs=50

    # fit the data
    clf = NadarayaWatson('gaussian', h=h_arr[i])
    clf.fit(z_sample[0:subs:, None], mu_sample[0:subs], dmu[0:subs])
    
    mu_sample_fit = clf.predict(z_sample[subs:, None])
    mu_fit = clf.predict(z[:, None])

    crossval1 = (np.sum((mu_sample_fit - mu_sample[subs:]) ** 2)
                / (len(mu_sample[subs:]) - 1))
    crossval[i]=crossval1

    ax.plot(z, mu_fit, '-', color='#DDDDDD')
    ax.plot(z, mu_true, '--', c='red')
    ax.errorbar(z_sample, mu_sample, dmu, fmt='.k', ecolor='gray', lw=1)

    #ax.text(0.5, 0.05, r"$\chi^2_{\rm dof} = %.2f$" % chi2_dof,
    #        ha='center', va='bottom', transform=ax.transAxes)

    ax.set_xlim(0.01, 1.8)
    ax.set_ylim(36.01, 48)
    ax.set_ylabel(r'$\mu$')
    ax.set_xlabel(r'$z$')


ax2.plot(h_arr, crossval, '-', color='#000000')
ax2.set_ylabel(r"$L_n(\^{f})$")
ax2.set_xlabel(r'$h$')

plt.show()
