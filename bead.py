# -*- coding: utf-8 -*-
"""
Created on Sat Sep 11 19:47:33 2021

@author: ojd34

The following set of functions takes a bead from DHPSFUAngleFast.m
and plots the localisations in 3D, fits a gaussian to find precision
in x,y,z.

"""

# Library imports and initialising.
import numpy as np
import matplotlib.pyplot as plot
import pandas as pd
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plot.rcParams.update({'font.size':10})
#from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
#mark_inset)
#import matplotlib.patches as patches
#from copy import copy
from matplotlib import cm
color = cm.rainbow(np.linspace(0, 1, 20))
from matplotlib.ticker import NullFormatter
nullfmt = NullFormatter()
import math


###########################################################################

# READ ME

# .3d files of SINGLE, CROPPED beads you wish to analyse.
# IMPORTANT: Use \\ to separate directories.

#datafile = r'C:\Users\ojd34\OneDrive - University of Cambridge\Desktop\analysis_2\50nm\numFrame_500_numPhotons_35000_radius_2.5e-08\DH.3d'

#xN yN zN iAll fAll data output into 3D file
#this code doesn't include photon intensities
###########################################################################
# File loading function.

def xy(datafile, fignum):
    data = np.loadtxt(datafile)
    xs = data[:, 0]
    ys = data[:, 1]
    zs = data[:, 2]
    iAll = data[:,3]
    fs = data[:, 4]
    def read_data(datafile):
        return [xs, ys, zs, fs]

    # if __name__ == '__main__':
    ###########################################################################
    # Load bead loc files (.3d) and set origin of scatterplots to zero.
    bead_one = read_data(datafile)
    bead_one[0] = bead_one[0] - np.mean(bead_one[0])
    bead_one[1] = bead_one[1] - np.mean(bead_one[1])

    ###########################################################################
    # Initialise settings for plots.
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.025
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.25]
    rect_histy = [left_h, bottom, 0.25, height]
    plot.figure(fignum, figsize=(4,4))
    axScatter = plot.axes(rect_scatter)
    axHistx = plot.axes(rect_histx)
    axHisty = plot.axes(rect_histy)

    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    ###########################################################################
    # Scatter plot of bead localisations (adjust x in
    # color[x] to change).

    colors= list(data[:, 4])
    axScatter.scatter(bead_one[0], bead_one[1], s=2, marker='o',
    c=colors, cmap="inferno", alpha=0.99)

    ###########################################################################
    # Set up extra plot settings and set up histograms of scatter plot positions.
    # axScatter.legend(loc=1, prop={'size': 8})
    axScatter.legend(loc=1, prop={'size': 8})
    
    binwidth = 2
    minx = np.min([bead_one[0]])
    maxx = np.max([bead_one[0]])
    miny = np.min([bead_one[1]])
    maxy = np.max([bead_one[1]])
    axScatter.set_xlim((minx, maxx))
    axScatter.set_ylim((miny, maxy))
    axScatter.set_xlabel('x (nm)')
    axScatter.set_ylabel('y (nm)')
    xbins = np.arange(minx, maxx + binwidth, binwidth)
    ybins = np.arange(miny, maxy + binwidth, binwidth)

    axHistx.hist(bead_one[0], bins=xbins, color=color[3], alpha=0.15)
    axHistx.set_xlim((minx, maxx))
    axHistx.set_ylabel('Frequency')
    axHistx.xaxis.set_ticks_position('none')
    axHisty.hist(bead_one[1], bins=ybins, color=color[3],
    orientation='horizontal', alpha=0.15)
    axHisty.set_ylim((miny, maxy))
    axHisty.set_xlabel('Frequency')
    axHisty.yaxis.set_ticks_position('none')
    ###########################################################################

    # Plot and fit the histograms of scatter plot positions to a Gaussian curve,
    #yielding the estimated bead/correction precision as sigma.
    from scipy.optimize import curve_fit
    def gauss(x, A, B, mu, sigma):
        return A + B * np.exp(-1.0 * (x - mu)**2 / (2 * sigma**2))
    x_hist, xbins_nn = np.histogram(bead_one[0], bins=xbins)
    y_hist, zbins_nn = np.histogram(bead_one[1], bins=ybins)
    xbinscenters = np.array([0.5 * (xbins[i] + xbins[i+1])
    for i in range(len(xbins)-1)])
    ybinscenters = np.array([0.5 * (ybins[i] + ybins[i+1])
    for i in range(len(ybins)-1)])
    ###########################################################################

    # Fit parameter outpust as poptx. The output is a list of the
    # estimated fit paramters for the paramters in the fitting function (gauss).
    # popt gives estimates for [A, B, mu, sigma] of the Gaussian function.
    # If the code doesn't fit your data, try changing the starting values of popt
    # in the lines below (make changes to p0).

    ############################### BEAD ONE ##################################
    #adjusts x_histogram (top)
    mean_x = sum(xbinscenters * x_hist) / sum(x_hist)
    sigma_x = np.sqrt(sum(x_hist * (xbinscenters - mean_x) **2) / sum(x_hist))
    poptx, pcovx = curve_fit(gauss, xdata=xbinscenters, ydata=x_hist, p0 =[min(x_hist), max(x_hist), mean_x, sigma_x])
    fitxs = np.linspace(minx, maxx, 1000)
    axHistx.plot(fitxs, gauss(fitxs, poptx[0], poptx[1], poptx[2], poptx[3]),
    linestyle='--', color='black', label='Gaussian Fit')

    #adjusts y_histogram(right)
    mean_y = sum(ybinscenters * y_hist) / sum(y_hist)
    sigma_y = np.sqrt(sum(y_hist * (ybinscenters - mean_y) **2) / sum(y_hist))
    popty, pcovy = curve_fit(gauss, xdata=ybinscenters, ydata=y_hist , p0 =[min(y_hist), max(y_hist), mean_y, sigma_y])
    fitys = np.linspace(miny, maxy, 1000)
    axHisty.plot(gauss(fitys, popty[0], popty[1], popty[2], popty[3]), fitys,
    linestyle='--', color='black', label='Gaussian Fit')

    ###########################################################################
    # Finish plots, save and print outputs for the fitting parameters.
    image_format = 'svg'                                # e.g .png, .svg, etc.
    resolution = 1200                                   # dpi
    image_name = datafile.replace(".","-") + '_xy.svg'  # the filename with dots replaced with dashes

    #plot.tight_layout()
    plot.savefig(image_name, format=image_format, dpi=resolution, bbox_inches='tight')
    plot.show()
    print('precision_x = ' + str(poptx[3]) + ' nm')
    print('precision_y = ' + str(popty[3]) + ' nm')

    ###########################################################################
    #Find FWHM of gaussian fits (pore measurements - superlocalise the pore)
    sigma_x = abs(poptx[3])
    sigma_y = abs(popty[3])
    confidence_x95 = 4*(abs(sigma_x))
    confidence_y95 = 4*(abs(sigma_y))
    print('95_confidence_x =', confidence_x95)
    print('95_confidence_y =', confidence_y95)
    confidence_x99 = 6*(abs(sigma_x))
    confidence_y99 = 6*(abs(sigma_y))
    print('99_confidence_x =', confidence_x99)
    print('99_confidence_y =', confidence_y99)
    print('\n')
    
    #FWHM calculations for truncating data x1 and x2 
    FWHM_x1 = (-(confidence_x99/2))
    # correctedFWHM_x1 = FWHM_x1 + np.mean(bead_one[0])
    # print('corrected_upper_bound_x =', correctedFWHM_x1)
    print('lowr_bound', FWHM_x1)
    
    FWHM_x2 =((confidence_x99/2))
    # correctedFWHM_x2 = FWHM_x2 + np.mean(bead_one[0])
    # print('corrected_lower_bound_x =', correctedFWHM_x2)
    print('upper_bound=', FWHM_x2)
 
    ###########################################################################
    # Save data (frame number, total photon number, precision in x, y and z; FWHM in x,y and z)
    maxfs = np.max(fs)
    meaniAll = np.mean(iAll)

    return maxfs, meaniAll, abs(poptx[3]), abs(popty[3]), confidence_x95, confidence_y95, confidence_x99, confidence_y99, FWHM_x1, FWHM_x2


def xz(datafile, fignum):
    data = np.loadtxt(datafile)
    xs = data[:, 0]
    ys = data[:, 2]
    zs = data[:, 2]
    iAll = data[:,3]
    fs = data[:, 4]
    def read_data(datafile):
        return [xs, ys, zs, fs]

    ###########################################################################
    # Load bead loc files (.3d) and set origin of scatterplots to zero.
    bead_one = read_data(datafile)
    bead_one[0] = bead_one[0] - np.mean(bead_one[0])
    bead_one[1] = bead_one[1] - np.mean(bead_one[1])

    ###########################################################################
    # Initialise settings for plots.
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.025
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.25]
    rect_histz = [left_h, bottom, 0.25, height]
    plot.figure(fignum, figsize=(4,4))
    axScatter = plot.axes(rect_scatter)
    axHistx = plot.axes(rect_histx)
    axHistz = plot.axes(rect_histz)

    axHistx.xaxis.set_major_formatter(nullfmt)
    axHistz.yaxis.set_major_formatter(nullfmt)

    ###########################################################################
    # Scatter plot of bead localisations (adjust x in
    # color[x] to change).

    colors= list(data[:, 4])
    axScatter.scatter(bead_one[0], bead_one[1], s=2, marker='o',
    c=colors, cmap="inferno", alpha=0.99)

    ###########################################################################
    # Set up extra plot settings and set up histograms of scatter plot positions.
    # axScatter.legend(loc=1, prop={'size': 8})
    axScatter.legend(loc=1, prop={'size': 8})
    binwidth = 2
    minx = np.min([bead_one[0]])
    maxx = np.max([bead_one[0]])
    minz = np.min([bead_one[1]])
    maxz = np.max([bead_one[1]])
    axScatter.set_xlim((minx, maxx))
    axScatter.set_ylim((minz, maxz))
    axScatter.set_xlabel('x (nm)')
    axScatter.set_ylabel('z (nm)')
    xbins = np.arange(minx, maxx + binwidth, binwidth)
    zbins = np.arange(minz, maxz + binwidth, binwidth)

    axHistx.hist(bead_one[0], bins=xbins, color=color[3], alpha=0.15)
    axHistx.set_xlim((minx, maxx))
    axHistx.set_ylabel('Frequency')
    axHistx.xaxis.set_ticks_position('none')
    axHistz.hist(bead_one[1], bins=zbins, color=color[3],
    orientation='horizontal', alpha=0.15)
    axHistz.set_ylim((minz, maxz))
    axHistz.set_xlabel('Frequency')
    axHistz.yaxis.set_ticks_position('none')
    ###########################################################################

    # Plot and fit the histograms of scatter plot positions to a Gaussian curve,
    #yielding the estimated bead/correction precision as sigma.
    from scipy.optimize import curve_fit
    def gauss(x, A, B, mu, sigma):
        return A + B * np.exp(-1.0 * (x - mu)**2 / (2 * sigma**2))
    x_hist, xbins_nn = np.histogram(bead_one[0], bins=xbins)
    z_hist, zbins_nn = np.histogram(bead_one[1], bins=zbins)
    xbinscenters = np.array([0.5 * (xbins[i] + xbins[i+1])
    for i in range(len(xbins)-1)])
    zbinscenters = np.array([0.5 * (zbins[i] + zbins[i+1])
    for i in range(len(zbins)-1)])
    ###########################################################################

    # Fit parameter outpust as poptx. The output is a list of the
    # estimated fit paramters for the paramters in the fitting function (gauss).
    # popt gives estimates for [A, B, mu, sigma] of the Gaussian function.
    # If the code doesn't fit your data, try changing the starting values of popt
    # in the lines below (make changes to p0).

    ############################### BEAD ONE ##################################
    #adjusts x_histogram (top)
    mean_x = sum(xbinscenters * x_hist) / sum(x_hist)
    sigma_x = np.sqrt(sum(x_hist * (xbinscenters - mean_x) **2) / sum(x_hist))
    poptx, pcovx = curve_fit(gauss, xdata=xbinscenters, ydata=x_hist, p0 =[min(x_hist), max(x_hist), mean_x, sigma_x])
    fitxs = np.linspace(minx, maxx, 1000)
    axHistx.plot(fitxs, gauss(fitxs, poptx[0], poptx[1], poptx[2], poptx[3]),
    linestyle='--', color='black', label='Gaussian Fit')

    #adjusts y_histogram( right)
    mean_z = sum(zbinscenters * z_hist) / sum(z_hist)
    sigma_z = np.sqrt(sum(z_hist * (zbinscenters - mean_z) **2) / sum(z_hist))
    poptz, pcovz = curve_fit(gauss, xdata=zbinscenters, ydata=z_hist, p0 =[min(z_hist), max(z_hist), mean_z, sigma_z])
    fitzs = np.linspace(minz, maxz, 1000)
    axHistz.plot(gauss(fitzs, poptz[0], poptz[1], poptz[2], poptz[3]), fitzs,
    linestyle='--', color='black', label='Gaussian Fit')

    ###########################################################################
    # Finish plots, save and print outputs for the fitting parameters.
    image_format = 'svg'                                # e.g .png, .svg, etc.
    resolution = 1200                                   # dpi
    image_name = datafile.replace(".","-") + '_xz.svg'  # the filename with dots replaced with dashes

    #plot.tight_layout()
    plot.savefig(image_name, format=image_format, dpi=resolution, bbox_inches='tight')
    plot.show()
    print('precision_x = ' + str(poptx[3]) + ' nm')
    print('precision_z = ' + str(poptz[3]) + ' nm')

    ###########################################################################
    #Find FWHM of gaussian fits (pore measurements - superlocalise the pore)
    sigma_x = abs(poptx[3])
    sigma_z = abs(poptz[3])
    
    confidence_x95 = 4*(abs(sigma_x))
    confidence_z95 = 4*(abs(sigma_z))
    print('95_confidence_x =', confidence_x95)
    print('95_confidence_z =', confidence_z95)
    confidence_x99 = 6*(abs(sigma_x))
    confidence_z99 = 6*(abs(sigma_z))
    print('99_confidence_x =', confidence_x99)
    print('99_confidence_z =', confidence_z99)
    print('\n')
    
    #FWHM calculations for truncating data x1 and x2 
    FWHM_z1 = (-(confidence_z99/2))
    print('lower_bound_z =', FWHM_z1)
    FWHM_z2 =(confidence_z99/2)   
    print('upper_bound_z =', FWHM_z2)

    ###########################################################################
    # Save data (frame number, total photon number, precision in x, y and z; FWHM in x,y and z)
    maxfs = np.max(fs)
    meaniAll = np.mean(iAll)

    return maxfs, meaniAll, abs(poptx[3]), abs(poptz[3]), confidence_x95, confidence_z95, confidence_x99, confidence_z99, FWHM_z1, FWHM_z2


def yz(datafile, fignum):
    ###########################################################################

    # READ ME

    # .3d files of SINGLE, CROPPED beads you wish to analyse.
    # IMPORTANT: Use \\ to separate directories.


    # datafile = 'C:\\Users\\ojd34\\OneDrive - University of Cambridge\\Desktop\\Precision DH-PSF\\5mW _3.results.3d'

    #xN yN zN iAll fAll data output into 3D file
    #this code doesn't include photon intensities
    ###########################################################################
    # File loading function.
    data = np.loadtxt(datafile)
    xs = data[:, 1]
    ys = data[:, 2]
    zs = data[:, 2]
    iAll = data[:,3]
    fs = data[:, 4]
    def read_data(datafile):
        return [xs, ys, zs, fs]

    ###########################################################################
    # Load bead loc files (.3d) and set origin of scatterplots to zero.
    bead_one = read_data(datafile)
    bead_one[0] = bead_one[0] - np.mean(bead_one[0])
    bead_one[1] = bead_one[1] - np.mean(bead_one[1])

    ###########################################################################
    # Initialise settings for plots.
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.025
    rect_scatter = [left, bottom, width, height]
    rect_histy = [left, bottom_h, width, 0.25]
    rect_histz = [left_h, bottom, 0.25, height]
    plot.figure(1, figsize=(4,4))
    axScatter = plot.axes(rect_scatter)
    axHisty = plot.axes(rect_histy)
    axHistz = plot.axes(rect_histz)

    axHisty.xaxis.set_major_formatter(nullfmt)
    axHistz.yaxis.set_major_formatter(nullfmt)

    ###########################################################################
        # Scatter plot of bead localisations (adjust x in
        # color[x] to change).

    colors= list(data[:, 4])
    axScatter.scatter(bead_one[0], bead_one[1], s=2, marker='o',
    c=colors, cmap="inferno", alpha=0.99)

    ###########################################################################
    # Set up extra plot settings and set up histograms of scatter plot positions.
    axScatter.legend(loc=1, prop={'size': 8})
    
    binwidth = 2
    miny = np.min([bead_one[0]])
    maxy = np.max([bead_one[0]])
    minz = np.min([bead_one[1]])
    maxz = np.max([bead_one[1]])
    axScatter.set_xlim((miny, maxy))
    axScatter.set_ylim((minz, maxz))
    axScatter.set_xlabel('y (nm)')
    axScatter.set_ylabel('z (nm)')
    ybins = np.arange(miny, maxy + binwidth, binwidth)
    zbins = np.arange(minz, maxz + binwidth, binwidth)

    axHisty.hist(bead_one[0], bins=ybins, color=color[3], alpha=0.15)
    axHisty.set_xlim((miny, maxy))
    axHisty.set_ylabel('Frequency')
    axHisty.xaxis.set_ticks_position('none')
    axHistz.hist(bead_one[1], bins=zbins, color=color[3],
    orientation='horizontal', alpha=0.15)
    axHistz.set_ylim((minz, maxz))
    axHistz.set_xlabel('Frequency')
    axHistz.yaxis.set_ticks_position('none')
    ###########################################################################

    # Plot and fit the histograms of scatter plot positions to a Gaussian curve,
    #yielding the estimated bead/correction precision as sigma.
    from scipy.optimize import curve_fit
    def gauss(x, A, B, mu, sigma):
        return A + B * np.exp(-1.0 * (x - mu)**2 / (2 * sigma**2))
    y_hist, ybins_nn = np.histogram(bead_one[0], bins=ybins)
    z_hist, zbins_nn = np.histogram(bead_one[1], bins=zbins)
    ybinscenters = np.array([0.5 * (ybins[i] + ybins[i+1])
    for i in range(len(ybins)-1)])
    zbinscenters = np.array([0.5 * (zbins[i] + zbins[i+1])
    for i in range(len(zbins)-1)])
    ###########################################################################

    # Fit parameter outpust as poptx. The output is a list of the
    # estimated fit paramters for the paramters in the fitting function (gauss).
    # popt gives estimates for [A, B, mu, sigma] of the Gaussian function.
    # If the code doesn't fit your data, try changing the starting values of popt
    # in the lines below (make changes to p0).

    ############################### BEAD ONE ##################################
    #adjusts y_histogram (top)
    mean_y = sum(ybinscenters * y_hist) / sum(y_hist)  
    sigma_y = np.sqrt(sum(y_hist * (ybinscenters - mean_y) **2) / sum(y_hist))
    popty, pcovy = curve_fit(gauss, xdata=ybinscenters, ydata=y_hist , p0 =[min(y_hist), max(y_hist), mean_y, sigma_y])
    fitys = np.linspace(miny, maxy, 1000)
    axHisty.plot(fitys, gauss(fitys, popty[0], popty[1], popty[2], popty[3]),
    linestyle='--', color='black')

    #adjusts z_histogram( right)
    mean_z = sum(zbinscenters * z_hist) / sum(z_hist)
    sigma_z = np.sqrt(sum(z_hist * (zbinscenters - mean_z) **2) / sum(z_hist))
    poptz, pcovz = curve_fit(gauss, xdata=zbinscenters, ydata=z_hist, p0 =[min(z_hist), max(z_hist), mean_z, sigma_z])
    fitzs = np.linspace(minz, maxz, 1000)
    axHistz.plot(gauss(fitzs, poptz[0], poptz[1], poptz[2], poptz[3]), fitzs,
    linestyle='--', color='black', label='Gaussian Fit')
    plot.legend()

    ###########################################################################
    # Finish plots, save plots and print outputs for the fitting parameters.
    image_format = 'svg'                                # e.g .png, .svg, etc.
    resolution = 1200                                   # dpi
    image_name = datafile.replace(".","-") + '_yz.svg'  # the filename with dots replaced with dashes

    #plot.tight_layout()
    plot.savefig(image_name, format=image_format, dpi=resolution, bbox_inches='tight')
    plot.show()
    print('precision_y = ' + str(popty[3]) + ' nm')
    print('precision_z = ' + str(poptz[3]) + ' nm')
    ###########################################################################
    #Find FWHM of gaussian fits (pore measurements - superlocalise the pore)
    sigma_y = abs(popty[3])
    sigma_z = abs(poptz[3])
    
    confidence_y95 = 4*(abs(sigma_y))
    confidence_z95 = 4*(abs(sigma_z))
    print('95_confidence_y =', confidence_y95)
    print('95_confidence_z =', confidence_z95)
    confidence_y99 = 6*(abs(sigma_y))
    confidence_z99 = 6*(abs(sigma_z))
    print('99_confidence_y =', confidence_y99)
    print('99_confidence_z =', confidence_z99)
    print('\n')
    
    #FWHM calculations for truncating data x1 and x2 
    FWHM_y1 = (-(confidence_y99/2))
    print('lower_bound_y =', FWHM_y1)
    FWHM_y2 =  ((confidence_y99/2))
    print('upper_bound_y =', FWHM_y2)

    # Save data (frame number, total photon number, precision in x, y and z; FWHM in x,y and z)
    maxfs = np.max(fs)
    meaniAll = np.mean(iAll)

    return maxfs, meaniAll, abs(popty[3]), abs(poptz[3]), confidence_y95, confidence_z95, confidence_y99, confidence_z99, FWHM_y1, FWHM_y2


















    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    