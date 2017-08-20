#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Frank Löffler
# license: CC BY 4.0 (http://creativecommons.org/licenses/by/4.0/)

from __future__ import print_function
import atexit
import numpy as np
import h5py, math, re, sys, os, os.path, fnmatch

if len(sys.argv) == 2 and sys.argv[1] == "--help":
  print("Usage: "+sys.argv[0]+" [filename]"+"""
Collect information about apparent horizon shape files within the current
directory and store them in hdf5 files. This assumes that the filename
scheme in AHFinderDirect (for ASCII output) has not been changed from the
default one.

The optional argument [filename] specifies another hdf5 file from which
timestep information will be read and used to create (if necessary dummy)
horizons for all timesteps in that file also in the apparent horizon output
files. The intend here is to allow visualization tools (like VisIt) to
automatically combine time sliders opon opening the specified file and the
created apparent horizon files.

The output files will be named AH_%d.[h5,xml], with %d being the number of each
horizon. The .h5 files contain the raw shape information, while the .xml files
describe its layout (and should be the one to open in a visualization tool).""")
  sys.exit(1)

h5alignfilename = None
if len(sys.argv) == 2:
  h5alignfilename = sys.argv[1]
  if not os.path.isfile(h5alignfilename):
    print("'"+h5alignfilename+"' is not a readable file.", file=sys.stderr)
    sys.exit(1)

# Class to store AH information into while it is read.
class AH(object):
  # Write stored xml to file - to be called when AH is destructed
  def write(self):
    self.xml += """
    </Grid>
  </Domain>
</Xdmf>"""
    with open(self.filename + '.xmf', 'w') as f:
      f.write(self.xml)

  # Constructor
  def __init__(self, filename):
    atexit.register(AH.write, self)
    self.filename   = filename
    self.h5filename = filename + '.h5'
    self.h5 = h5py.File(self.h5filename, 'w')
    self.xml = """<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0">
  <Domain>
    <Grid Name="AHs" GridType="Collection" CollectionType="Temporal">"""

  # Read AH data from file 'filename'.
  def get_ah(self, filename):
    # Get shape information from file
    N_rho   = None
    N_sigma = None
    with open(filename, 'r') as f:
      for line in f:
        if line[0:10] == '# N_rho = ':
          regex = re.match(r'# N_rho = (\d+)', line)
          if (regex):
            N_rho = int(regex.group(1))
        if line[0:12] == '# N_sigma = ':
          regex = re.match(r'# N_sigma = (\d+)', line)
          if (regex):
            N_sigma = int(regex.group(1))
        if (N_rho != None and N_sigma != None):
          break
    # Bail out if this doesn't look like a file we expect
    if not N_rho or not N_sigma:
      print("Could not find N_rho and N_sigma. Is this a uncorrupted AHFinderdirect file?",
            file=sys.stderr)
      sys.exit(1)
    if N_rho != N_sigma:
      print("N_rho (%d) != N_sigma (%d). This is completely untested and probably does not "
            "work" % (N_rho, N_sigma), file=sys.stderr)
      sys.exit(1)
    # read actual data
    (x,y,z) = np.genfromtxt(filename, usecols=(3,4,5), comments='#', unpack=True)
    x = x.reshape((6, N_rho, N_sigma))
    y = y.reshape((6, N_rho, N_sigma))
    z = z.reshape((6, N_rho, N_sigma))
    return (N_rho, N_sigma, x, y, z)

  # Append data from AH in filename (iteration it) to xdmf
  def append(self, filename, it):
    # If we don't have an AH for this timestep, add a dummy mesh.
    if (filename == None):
      (N_rho, N_sigma, x, y, z) = (1, 1, np.zeros((6,1,1)),
                                         np.zeros((6,1,1)),
                                         np.zeros((6,1,1)))
    else:
      (N_rho, N_sigma, x, y, z) = self.get_ah(filename)
    self.xml += """
      <Grid Name="AH" GridType="Collection">
        <Time Value="%(it)d" />""" % {'it': it}
    # copy data into hdf5 file and save corresponding xml
    for i in range(6):
      group = '/it_%d/h_%d/p_%d/' % (it, 1, i+1)
      self.h5[group + 'x'] = x[i,:,:]
      self.h5[group + 'y'] = y[i,:,:]
      self.h5[group + 'z'] = z[i,:,:]
      hlink = '/h_%d/it_%d/p_%d/' % (1, it, i+1)
      self.h5[hlink + 'x'] = h5py.SoftLink(group + 'x')
      self.h5[hlink + 'y'] = h5py.SoftLink(group + 'y')
      self.h5[hlink + 'z'] = h5py.SoftLink(group + 'z')

      self.xml += """
        <Grid Name="AH_mesh%(i)d" GridType="Uniform">
          <Topology TopologyType="3DSMesh" NumberOfElements="%(xy)s"/>
          <Geometry GeometryType="X_Y_Z">
            <DataItem Dimensions="%(xy)s" NumberType="Float" Precision="8" Format="HDF">
              %(h5fn)s:%(group)sx
            </DataItem>
            <DataItem Dimensions="%(xy)s" NumberType="Float" Precision="8" Format="HDF">
              %(h5fn)s:%(group)sy
            </DataItem>
            <DataItem Dimensions="%(xy)s" NumberType="Float" Precision="8" Format="HDF">
              %(h5fn)s:%(group)sz
            </DataItem>
          </Geometry>
        </Grid>""" % {'i': i, 'i1': i+1, 'xy': "%d %d" % (N_rho, N_sigma),
                      'h5fn': self.h5filename, 'group': group}
    self.xml += """
      </Grid>"""

# figure out which iterations are present in related hdf5 output
h5_align_its = {}
if (h5alignfilename != None):
  with h5py.File(h5alignfilename, 'r') as h5_align:
    for item in h5_align.keys():
      reg_it = re.search(r'\S+ it=(\d+) ', item)
      if reg_it:
        h5_align_its[int(reg_it.group(1))] = 1
h5_align_its = np.asarray(sorted(h5_align_its.keys()))

# find number of horizons
hs = {}
for root, subdirs, files in os.walk('.'):
  for filename in files:
    if fnmatch.fnmatch(filename, 'h.t*.ah*.gp'):
      h = int(re.search('h\.t\d+\.ah(\d+)\.gp', filename).group(1))
      hs[h] = 1
hs = np.asarray(sorted(hs.keys()))

# find timesteps for each horizon
ah_align_fns = {}
ah_align_its = {}
for h in hs:
  ah_align_fns[h] = {}
  for root, subdirs, files in os.walk('.'):
    for filename in files:
      if fnmatch.fnmatch(filename, 'h.t*.ah%d.gp' % (h)):
        it = int(re.search('h\.t(\d+)\.ah%d\.gp' % (h), filename).group(1))
        ah_align_fns[h][it] = os.path.join(root, filename)
  ah_align_its[h] = np.array(sorted(ah_align_fns[h].keys()))

# Generation of intersecting iterations, and some info output
if (h5alignfilename != None):
  print("Found shapes for %d horizon(s), and %d iterations in hdf5 output." %
        (hs.size, h5_align_its.size))
else:
  print("Found shapes for %d horizon(s)." %
        (hs.size))
h5_align_interect = {}
for h in hs:
  h5_align_interect[h] = np.intersect1d(h5_align_its, ah_align_its[h], assume_unique=True)
  if (h5alignfilename != None):
    print("horizon %d: %d timesteps, out of which %d intersect with hdf5." %
          (h, ah_align_its[h].size, h5_align_interect[h].size))
  else:
    print("horizon %d: %d timesteps." % (h, ah_align_its[h].size))

if (h5alignfilename != None):
  # Add horizon for each iteration in h5 file: write dummy if not present
  for h in hs:
    ah = AH("AH_%d" % (h))
    for it in h5_align_its:
      if it in h5_align_interect[h]:
        ah.append(ah_align_fns[h][it], it)
      else:
        ah.append(None, it)
else:
  # add all iterations to h5 file
  for h in hs:
    ah = AH("AH_%d" % (h))
    for it in ah_align_fns[h].keys():
      ah.append(ah_align_fns[h][it], it)

