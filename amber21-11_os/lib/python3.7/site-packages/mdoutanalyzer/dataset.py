#!/usr/bin/env python

# Make sure all division is floating division


import numpy as np
import math
import sys
   
class DatasetError(Exception):
   """ If you gave a bad parameter... """

class DataSet(np.ndarray):
   """ Derived class from Numpy array class """

   def __array_finalize__(self, obj):
      """ Set up a counter so we can add to this """
      self.counter = 0
#     np.ndarray.__array_finalize__(self, obj)

   def KullbackLeibler(self, interval, ofile=sys.stdout):
      r"""
      Performs Kullback-Liebler analysis on the data set by histogramming
      portions of the data.  New histograms are collected every 'interval'
      steps.

      Kullback-Leibler divergence compares two distributions (in this case
      the one at time t and the final distribution) according to the formula
                i=N  /        P(i) \
      D(P||Q) = sum | P(i) ln ----  |
                i=1  \        Q(i) /
      """
      close_file_after = False
      if type(ofile).__name__ == 'str':
         close_file_after = True
         try:
            ofile = open(ofile, 'w')
         except IOError:
            raise DatasetError("Could not open %s for writing!" % ofile)

      if not hasattr(ofile, 'write'):
         raise TypeError("ofile must be of type 'str' or derived from 'file'!")

      if interval > self.size:
         raise DatasetError("interval (%d) is larger than my size (%d)!" %
                            interval, self.size)
      
      # Make sure we set up our histogram preferences
      if not hasattr(self, 'hmin'):
         raise DatasetError("Histogram parameters are not set!")

      # Calculate our total histogram if we haven't already
      if not hasattr(self, 'tot_dist'):
         self._finalhist()

      if type(interval).__name__ != 'int':
         raise TypeError("interval expected to be an integer!")

      for i in range(1, (self.size-1) // interval+1):
         hist = np.histogram(self[:i*interval],
                             bins=self.nbins,
                             range=(self.hmin, self.hmax),
                             weights=None,
                             density=self.norm)
         ofile.write("%d %g\n" % (i, _kull_leib(hist[0], self.tot_dist[0])))
      hist = np.histogram(self[:],
                          bins=self.nbins,
                          range=(self.hmin, self.hmax),
                          weights=None,
                          density=self.norm)
      ofile.write("%d %g\n" % (i+1, _kull_leib(hist[0], self.tot_dist[0])))

      if close_file_after: ofile.close()

   def set_hist_params(self, hmin=None, hmax=None, nbins=None, spacing=None,
                       norm=True):
      """ 
      Set the histogram parameters -- convert everything to max/min and nbins
      """

      # hmin/hmax are set from min/max unless explicitly specified
      if hmin is None:
         self.hmin = self.min()
      else:
         self.hmin = hmin
      if hmax is None:
         self.hmax = self.max()
      else:
         self.hmax = hmax

      # nbins and spacing are exclusive
      if nbins is not None and spacing is not None:
         raise DatasetError("Cannot specify both nbins and spacing!")

      # Scott's Choice for spacing if neither are specified
      if spacing is None:
         spacing = 3.5 * self.std() / self.size ** (1/3)

      # Take our bins if specified, or use spacing (which has a default of
      # Scott's Choice set above if it's not already set)
      if nbins is not None:
         self.nbins = nbins
      else:
         self.nbins = int(math.ceil((self.hmax - self.hmin) / spacing))

      # Do we want to normalize this histogram?
      self.norm = norm

   def _finalhist(self):
      """ This analyzes the 1dhistogram of the full data set """
      self.tot_dist = np.histogram(self, 
                                   bins=self.nbins,
                                   range=(self.hmin, self.hmax),
                                   weights=None,
                                   density=self.norm)
   def append(self, val):
      """ Adds a new element to the array """
      # Assume it is a float, to save time
      self.resize(self.size+1, refcheck=False)
      self[self.size-1] = val

   def add_value(self, val):
      """ Add a float to my array """
      if not isinstance(val, float):
         raise TypeError('add_value expects a float!')
      if self.counter >= self.shape[0]:
         raise IndexError('Cannot add_value to DataSet, out of bounds (%d)!' % 
                          self.counter)
      self[self.counter] = val
      self.counter += 1

   def truncate(self):
      """ Get rid of everything after counter """
      return np.resize(self, self.counter)

def _kull_leib(hist1, hist2):
   """ 
   Integrates two histograms hist1 and hist2 according to the formula shown
   above for Kullback-Leibler, where hist1 is P and hist2 is Q
   """
   if hist1.size != hist2.size:
      raise DatasetError("Cannot integrate distributions of different sizes!")

   runsum = 0.0
   for i in range(hist1.size):
      if hist1[i] == 0 or hist2[i] == 0: continue
      runsum += hist1[i] * math.log(hist1[i]/hist2[i])

   return runsum

def load_from_file(infile, column):
   """ Loads a DataSet from a file """
   ret = DataSet(0)
   for line in infile:
      try:
         val = float(line.split()[column-1])
      except IndexError: continue
      except ValueError: continue
      ret.resize(ret.size+1, refcheck=False)
      ret[ret.size-1] = val

   return ret
