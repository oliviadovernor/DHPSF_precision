# -*- coding: utf-8 -*-
"""
Created on Sat Sep 11 19:45:47 2021

@author: ojd34
"""

##############################################################################
# Importing libraries

import numpy as np
import os
import csv
import bead
import histogram_FWHM
import xlsxwriter
from scipy.optimize._lsq.common import print_header_linear

##############################################################################
# Setting root directory - should be folder containing main.py and data to be
# analysed in adjacent folder

rootdir = r'C:\Users\ojd34\OneDrive - University of Cambridge\Desktop\Precision data'

##############################################################################

# truncating data
  #  1. open up 3d file 
  #  2. convert to csv 
  #  3. take into account limits upper and lower bounds - delete data outside bounds 
  #  4. save the csv to same file _ truncated name 

def truncate(datafile, xLB, xUB, zLB, zUB, yLB, yUB):
# for subdir, dirs, files in os.walk(rootdir):
#     # find .3d files
#     for file in files:
#         if file.endswith(".3d"):
#             filepath = os.path.join(subdir, file) # copy path
            
    data = np.loadtxt(datafile)
    xs = data[:, 0]
    ys = data[:, 1]
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

    #Undoing 0,0 origin setting of 2D plots (simulated data)
    
    mean_oldx = np.mean(xs)
    #originalx = mean_oldx + bead_one[0]
    
    mean_oldy = np.mean(ys)
    #originaly = mean_oldy + bead_one[1]
    
    mean_oldz = np.mean(zs)
    #originalz = mean_oldz + bead_one[2]

 #1.  pick from list of new values (y set gregs set or bead 1) which are inside the limits worked out for upper and lower bounds FWHM
 
#    print(len(xs))
    
    count = 0
    for row in range(len(data)):
         if ((xLB+mean_oldx <= bead_one[0][row]+mean_oldx <= xUB+mean_oldx) or (yLB+mean_oldy <= bead_one[1][row]+mean_oldy <= yUB+mean_oldy) or (zLB+mean_oldz <= bead_one[2][row]+mean_oldz <= zUB+mean_oldz) or ((xLB+mean_oldx <= bead_one[0][row]+mean_oldx <= xUB+mean_oldx) and (yLB+mean_oldy <= bead_one[1][row]+mean_oldy <= yUB+mean_oldy) and (zLB+mean_oldz <= bead_one[2][row]+mean_oldz <= zUB+mean_oldz))):
             print(str(xLB+mean_oldx) + " <= " + str(bead_one[0][row]+mean_oldx) + " <= " + str(xUB+mean_oldx))
             print(str(xLB) + " <= " + str(bead_one[0][row]) + " <= " + str(xUB))
             count += 1
             
             truncDataFilename = file.replace(".","-") + '_trunc_FWHMdata.csv'  # the filename with dots replaced with dashes
             truncDataFilepath = os.path.join(subdir, truncDataFilename) # create path
            
             with open(truncDataFilepath, 'a', newline='') as newData:
                 writeNewData = csv.writer(newData, delimiter = ',', quotechar = '"', quoting = csv.QUOTE_MINIMAL)
                 writeNewData.writerow([xs[row], ys[row], zs[row], iAll[row], fs[row]])
    
    print(count)

##############################################################################
# Initialising .xlsx output file by creating it and inserting col headers

# Workbook() takes one, non-optional, argument
# which is the filename that we want to create.
datalog = xlsxwriter.Workbook('log.xlsx')

# The workbook object is then used to add new worksheet
sheet1 = datalog.add_worksheet()

# Use the worksheet object to write data
sheet1.write('A1', 'Number')    # Column title
sheet1.write('B1', 'Filename')  # validation column - True/1 if frame numbers match, False/0 otherwise

sheet1.write('C1', 'Total Frame number')    # Column title
sheet1.write('D1', 'Valid')                 # validation column - True/1 if frame numbers match, False/0 otherwise

sheet1.write('E1', 'Average Total Photon number')   # Column title
sheet1.write('F1', 'Valid')                         # validation column - True/1 if photon numbers match, False/0 otherwise

sheet1.write('G1', 'precision_x')   # Column title
sheet1.write('H1', 'Valid')         # validation column - True/1 if all x precisions match, False/0 otherwise

sheet1.write('I1', 'precision_y')   # Column title
sheet1.write('J1', 'Valid')         # validation column - True/1 if all y precisions match, False/0 otherwise

sheet1.write('K1', 'precision_z')   # Column title
sheet1.write('L1', 'Valid')         # validation column - True/1 if all z precisions match, False/0 otherwise

sheet1.write('M1', '95_confidence_x')    # Column title
sheet1.write('N1', 'Valid')     # validation column - True/1 if 95% in x match, False/0 otherwise

sheet1.write('O1', '95_confidence_y')    # Column title
sheet1.write('P1', 'Valid')     # validation column - True/1 if 95% in y match, False/0 otherwise

sheet1.write('Q1', '95_confidence_z')    # Column title
sheet1.write('R1', 'Valid')     # validation column - True/1 if 95% in z match, False/0 otherwise

sheet1.write('S1', '99_confidence_x')    # Column title
sheet1.write('T1', 'Valid')     # validation column - True/1 if 99% in x match, False/0 otherwise

sheet1.write('U1', '99_confidence_y')    # Column title
sheet1.write('V1', 'Valid')     # validation column - True/1 if 99% in y match, False/0 otherwise

sheet1.write('W1', '99_confidence_z')    # Column title
sheet1.write('X1', 'Valid')     # validation column - True/1 if 99% in z match, False/0 otherwise

##############################################################################
# Iterating through .3d files

filecount = 0 # initialising the file counter/row number in log.xlsx
plotcount = 0 # initialising the plot index/figure number in bead.py

# search folders in workspace
for subdir, dirs, files in os.walk(rootdir):
    # find .3d files
    for file in files:
        if file.endswith(".3d"):
            filepath = os.path.join(subdir, file) # copy path
            print(filepath)

            filecount +=1   # filecount now actually represents the number of .3d files and is used to index the rows of the xlsx
            print('\n')

            plotcount += 1  # indexes the plots
            maxfs_xy, meaniAll_xy, xPrscn_xy, yPrscn_xy, confidence_x95_xy, confidence_y95_xy, confidence_x99_xy, confidence_y99_xy, xyFWHM_x1, xyFWHM_x2 = bead.xy(filepath,plotcount) # unpack return parameters of bead.xy()

            plotcount += 1  # indexes the plots
            maxfs_xz, meaniAll_xz, xPrscn_xz, zPrscn_xz, confidence_x95_xz, confidence_z95_xz, confidence_x99_xz, confidence_z99_xz, xzFWHM_z1, xzFWHM_z2 = bead.xz(filepath,plotcount) # unpack return parameters of bead.xz()

            plotcount += 1  # indexes the plots
            maxfs_yz, meaniAll_yz, yPrscn_yz, zPrscn_yz, confidence_y95_yz, confidence_z95_yz, confidence_y99_yz, confidence_z99_yz, yzFWHM_y1, yzFWHM_y2 = bead.yz(filepath,plotcount) # unpack return parameters of bead.yz()

            # construct the row of data
            row = [ filecount, file,
                        maxfs_xy, maxfs_xy == maxfs_xz == maxfs_yz,             # [maxfs,       check maxfs matches,
                        meaniAll_xy, meaniAll_xy == meaniAll_xy == meaniAll_xy, #  meaniAll,    check meaniAll matches
                        xPrscn_xy, xPrscn_xy == xPrscn_xz,                      #  x precision, check x precision matches
                        yPrscn_xy, yPrscn_xy == yPrscn_yz,                      #  y precision, check y precision matches
                        zPrscn_xz, zPrscn_xz == zPrscn_yz,                      #  z precision, check z precision matches
                        confidence_x95_xy, confidence_x95_xy == confidence_x95_xz,                         #  x FWHM,      check x FWHM matches
                        confidence_y95_xy, confidence_y95_xy == confidence_y95_yz,                         #  y FWHM,      check y FWHM matches
                        confidence_z95_xz, confidence_z95_xz == confidence_z95_yz, 
                        confidence_x99_xy, confidence_x99_xy == confidence_x99_xz,                         #  x FWHM,      check x FWHM matches
                        confidence_y99_xy, confidence_y99_xy == confidence_y99_yz,                         #  y FWHM,      check y FWHM matches
                        confidence_z99_xz, confidence_z99_xz == confidence_z99_yz]                         #  z FWHM,      check z FWHM matches]

            # insert the row of data into the .xlsx file
            for column, item in enumerate(row):
                    sheet1.write(filecount,column,item)
                    
            #print("Truncating...")
            #truncate(filepath, xyFWHM_x1, xyFWHM_x2, xzFWHM_z1, xzFWHM_z2, yzFWHM_y1, yzFWHM_y2)

# Close the Excel file
datalog.close()

# Mean and distribution of pore sizes: 
# Calls function which plots x, y and z FWHM data into a histogram 
# mean and st dev data collected for each (mean = average pore size in each dimension, st dev = distribution info)
histogram_FWHM.x_confidence()
histogram_FWHM.y_confidence()
histogram_FWHM.z_confidence()

# x, y and z (lateral and axial) precision line plot:
    
            


