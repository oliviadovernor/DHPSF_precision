# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 22:15:22 2021

@author: ojd34
"""
#Plot histogram pf FWHMs


import pandas as pd
from matplotlib import pyplot as plt
import statistics

def x_confidence ():

    # Read spreadsheet created in main containing pore data in x, y and z
    df = pd.read_excel('log.xlsx', sheet_name = 'Sheet1', engine = 'openpyxl')
    data_points = df['95_confidence_x'].tolist()
    
    # Adjust bins/sections shown within the histogram 
    n_bins =20
    
    # Plot histogram
    fig, ax = plt.subplots()
    arr = ax.hist(data_points, ec='black', bins=n_bins)
    
    # Adjust x axis range and section (bin) size
    plt.xticks(range(0, 100, n_bins))
    
    # Add labels to indicate how many datapoints are in subcolumn    
    labels = [arr[0][i] for i in range(len(arr[0]))]
    
    # Set y limit with +1 height to top (white space)
    ax.set_ylim(0, max(labels) + 5)
    
    rects = ax.patches
    
    for rect, label in zip(rects, labels): 
        height = rect.get_height()
        ax.text(
            rect.get_x() + rect.get_width() /2, 
            height + 0.01,
            int(label), 
            ha='center', 
            va='bottom'
            )
    
    # Calculate mean and standard deviation (sd and 2*sd) for x data 
    m_x = statistics.mean(data_points)
    sd_x = statistics.stdev(data_points)
    
    plt.axvline(m_x, color='k', linestyle ='dashed')
    plt.axvline(m_x + sd_x, color='y', linestyle ='dashed')
    plt.axvline(m_x - sd_x, color='y', linestyle ='dashed')
    plt.axvline(m_x + 2*sd_x, color='y', linestyle ='dashed')
    plt.axvline(m_x - 2*sd_x, color='y', linestyle ='dashed')
    
    plt.title('95% Confidence histogram_x')
    plt.show()
    print('mean 95% confidence_x =', m_x)
    print('stand deviation of 95% results =', sd_x)
    
    return m_x ,sd_x
  

def y_confidence ():

    # Read spreadsheet created in main containing pore data in x, y and z
    df = pd.read_excel('log.xlsx', sheet_name = 'Sheet1', engine = 'openpyxl')
    data_points = df['95_confidence_y'].tolist()
    
    # Adjust bins/sections shown within the histogram
    n_bins = 20
    
    # Plot histogram
    fig, ax = plt.subplots()
    arr = ax.hist(data_points, ec='black', bins=n_bins)
    
    # Adjust x axis range and section (bin) size
    plt.xticks(range(0, 120, n_bins))
    
    # Add labels to indicate how many datapoints are in subcolumn    
    labels = [arr[0][i] for i in range(len(arr[0]))]
    
    # Set y limit with +1 height to top (white space)
    ax.set_ylim(0, max(labels) + 5)
    
    rects = ax.patches
    
    for rect, label in zip(rects, labels): 
        height = rect.get_height()
        ax.text(
            rect.get_x() + rect.get_width() /2, 
            height + 0.01,
            int(label), 
            ha='center', 
            va='bottom'
            )
    
    # Calculate mean and standard deviation (sd and 2*sd) for x data 
    m_y = statistics.mean(data_points)
    sd_y = statistics.stdev(data_points)
    
    plt.axvline(m_y, color='k', linestyle ='dashed')
    plt.axvline(m_y + sd_y, color='g', linestyle ='dashed')
    plt.axvline(m_y - sd_y, color='g', linestyle ='dashed')
    plt.axvline(m_y + 2*sd_y, color='g', linestyle ='dashed')
    plt.axvline(m_y - 2*sd_y, color='g', linestyle ='dashed')
    
    plt.title('95% Confidence histogram_y')
    plt.show()
    print('mean 95% confidence_y =', m_y)
    print('stand deviation of 95% results =', sd_y)
    
    return m_y ,sd_y

def z_confidence ():

    # Read spreadsheet created in main containing pore data in x, y and z
    df = pd.read_excel('log.xlsx', sheet_name = 'Sheet1', engine = 'openpyxl')
    data_points = df['95_confidence_z'].tolist()
    
    # Adjust bins/sections shown within the histogram
    n_bins = 20
    
    # Plot histogram
    fig, ax = plt.subplots()
    arr = ax.hist(data_points, ec='black', bins=n_bins)
    
    # Adjust x axis range and section (bin) size
    plt.xticks(range(0, 250, n_bins))
    
    # Add labels to indicate how many datapoints are in subcolumn
    labels = [arr[0][i] for i in range(len(arr[0]))]
    
    # Set y limit with +1 height to top (white space)
    ax.set_ylim(0, max(labels) + 5)
    
    
    rects = ax.patches
    
    for rect, label in zip(rects, labels): 
        height = rect.get_height()
        ax.text(
            rect.get_x() + rect.get_width() /2, 
            height + 0.01,
            int(label), 
            ha='center', 
            va='bottom'
            )
    
    # Calculate mean and standard deviation (sd and 2*sd) for x data 
    m_z = statistics.mean(data_points)
    sd_z = statistics.stdev(data_points)
    
    plt.axvline(m_z, color='k', linestyle ='dashed')
    plt.axvline(m_z + sd_z, color='r', linestyle ='dashed')
    plt.axvline(m_z - sd_z, color='r', linestyle ='dashed')
    plt.axvline(m_z + 2*sd_z, color='r', linestyle ='dashed')
    plt.axvline(m_z - 2*sd_z, color='r', linestyle ='dashed')
    
    plt.title('95% Confidence histogram_z')
    plt.show()
    print('mean 95% confidence_z =', m_z)
    print('stand deviation of 95% results =', sd_z)
    
    return m_z ,sd_z
        

    
