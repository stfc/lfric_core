#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################
'''
Basic python script to plot the x-z profile minus a constant state of 300
from a Dynamo output file.

This version takes nodal format output files and
interpolates onto a regular grid.

Levels are determined from the data

This version stitches together a directory of files
and extracts all levels so it can work in the serial
case where there is one file or the parallel case where
there is a file for each processor.

'''

import numpy as np
# Need to set a non-interactive backend for suites
import matplotlib
matplotlib.use('Agg')

import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from scipy.interpolate import griddata

import math

import glob
import sys


# Make an empty list to hold the levels we find in the data
levels = []

# Set up some empty lists for
# x,y coordinates and value
x = []
y = []
z = []

def process_file_list(filestem):
  # get the list of files to stitch together
  dirlist = glob.glob(filestem)
 
  for f in dirlist:

    print "processing file ", f
    fo = open(f, "r")

    # Step through all lines in the file, split the lines
    # and where the level matches the specifed one, append 
    # data to appropriate list 
    for strline in fo:
       strsplit = strline.split()
       # check we got a valid data line
       if (len(strsplit) == 5):
          # get the level
          level = float(strsplit[3])
          # Is the level already in the levels list?
          if (level in levels):
             # If it is then append the data into the correct list
             x[levels.index(level)].append(float(strsplit[0]))
             y[levels.index(level)].append(float(strsplit[1]))
             z[levels.index(level)].append(float(strsplit[4]))
          else:
             # add the level to the levels list and append
             # corresponding empty lists to x, y and z lists
             levels.append(level)
             x.append([])
             y.append([])
             z.append([])
             # ...and then append the data
             x[levels.index(level)].append(float(strsplit[0]))
             y[levels.index(level)].append(float(strsplit[1]))
             z[levels.index(level)].append(float(strsplit[4]))

    fo.close()
       
def make_figure(plotpath, field, timestep):

  slice_fig = plt.figure(figsize=(15,10))
  # get min and max of x,y data for plot axes
  xmin =  min(x[0])
  xmax = max(x[0])
  ymin =  min(y[0])
  ymax = max(y[0])
  zmin = 0.0
  zmax = 6400.0

  r2d = 1.0/1000.0;
  ny = 300
  nx = 2
  nz = len(levels)

  #create 2D plot
  x2d = 0.0
  z2d = np.linspace(zmin, zmax, nz)
  y2d = np.linspace(ymin, ymax, ny)
  xi, yi = np.meshgrid(x2d, y2d)    
  zi = np.zeros([ny,1,len(levels)])
  for p in xrange(len(levels)):
    zi[:,:,p] = griddata((np.asarray(x[p]), np.asarray(y[p])), np.asarray(z[p]), (xi, yi), method='linear') 

  yi, xi = np.meshgrid(z2d, y2d) 
  dz = np.zeros([ny,len(levels)])
  for i in range(ny):
    dz[i,:] = zi[i,0,:] - 300.0


  cc = np.linspace(-16,-1,16)
  cf = plt.contourf(xi *r2d, yi * r2d, dz, cc)
  cl = plt.contour(xi * r2d, yi*r2d, dz, cc, linewidths=0.5,colors='k')
  plt.axis([0, 16, 0, 5])
  plt.title('max: %.6s, min: %.6s'%(np.max(dz),np.min(dz)))

  out_file_name = plotpath + "/" + 'straka_y' + "_" + timestep +  ".png"
  slice_fig.savefig(out_file_name , bbox_inches='tight')

if __name__ == "__main__":
  
  try:
    datapath, fields, timesteps, plotpath = sys.argv[1:5]
  except ValueError:
    print("Usage: {0} <datapath> <field_names> <timestep_list> <plotpath>".format(sys.argv[0]))
    exit(1)

  # Split out the list of fields
  field_list = fields.split(':')

  # Split out the list of timesteps
  ts_list = timesteps.split(':')

  for field in field_list:

    for ts in ts_list:

      # clear the lists in between plots
      del levels[:]
      del x[:]
      del y[:]
      del z[:]

      filestem =  datapath + "/diagDynamo_nodal_" + field + "_" + ts + "*"

      process_file_list(filestem)
      # Only try to plot if we found some files for this timestep
      if len(levels) > 0:
        make_figure(plotpath,field, ts)

