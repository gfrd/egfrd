#!/usr/env python

import math
import numpy

from single import *
from pair import *
from multi import *

__all__ = [ 'Histogram3D', 'DomainsHistogramCollection' ]


class Histogram3D(object):

    def __init__(self, nbins, dimensions, name='noname'):

        self.name       = str(name)

        self.nbins      = nbins
        self.dimensions = dimensions
        self.histogram  = numpy.zeros(nbins)
        self.total_cnt  = 0

        self.binsizes   = self.calculate_binsizes()

        self.normalized = False


    def calculate_binsizes(self):

        binsizes = numpy.zeros(3)

        for d in [0, 1, 2]:

            binsizes[d] = 1.0 * self.dimensions[d] / self.nbins[d]

        if not all([bs>0 for bs in binsizes]):
            raise HistogramError('At least one bin size was zero in histogram %s.' % self.__str__() )

        return binsizes


    def bin(self, data):

        index = [-1, -1, -1]

        for d in [0, 1, 2]:

            b = int(data[d] / self.binsizes[d])

            if b >= 0 and b < self.nbins[d]:

                    index[d] = b

        ilist = [index[0], index[1], index[2]]

        if all([i>=0 and i<=self.nbins for i in ilist]):
            
            [ix, iy, iz] = ilist
            self.histogram[ix][iy][iz] = self.histogram[ix][iy][iz] + 1

            self.total_cnt = self.total_cnt + 1


    def normalize(self):

        if not all([bs>0 for bs in self.binsizes]):
            raise HistogramError('At least one bin size was zero in histogram %s.' % self.__str__() )

        if not self.total_cnt > 0:
            raise HistogramError('No counts in histogram %s, cannot normalize!' % self.__str__() )
        
        binvolume  = self.binsizes[0] * self.binsizes[1] * self.binsizes[2]
        normfactor = self.total_cnt * binvolume

        for i0 in range(0, self.nbins[0]):
          for i1 in range(0, self.nbins[1]):
            for i2 in range(0, self.nbins[2]):

                self.histogram[i0][i1][i2] = self.histogram[i0][i1][i2] / normfactor


    def integrate(self):
        
        binvolume  = self.binsizes[0] * self.binsizes[1] * self.binsizes[2]

        I = 0.0
        for i0 in range(0, self.nbins[0]):
          for i1 in range(0, self.nbins[1]):
            for i2 in range(0, self.nbins[2]):

                I = I + self.histogram[i0][i1][i2] * binvolume

        return I


    def check(self):

        I = self.integrate()

        if not feq(I, 1.0):

            raise HistogramError( 'Integrated density does not equal 1 in histogram %s; integral = %s.' % (self.__str__(), I) )


    def __str__(self):

        if self.name != '':

            return self.name

        else:
            return 'Histogram3D'

    

class DomainsHistogramCollection(object):

    def __init__(self, world, nbins, name='noname'):

        self.name = name

        self.world = world
        self.ws = world.world_size

        self.nbins = nbins

        self.SinglesHistogram   = Histogram3D([self.nbins, self.nbins, self.nbins], [self.ws, self.ws, self.ws], name='SinglesHistogram')
        self.PairsHistogram     = Histogram3D([self.nbins, self.nbins, self.nbins], [self.ws, self.ws, self.ws], name='PairsHistogram')
        self.MultisHistogram    = Histogram3D([self.nbins, self.nbins, self.nbins], [self.ws, self.ws, self.ws], name='MultisHistogram')
        self.NonMultisHistogram = Histogram3D([self.nbins, self.nbins, self.nbins], [self.ws, self.ws, self.ws], name='NonMultisHistogram')


    def bin_domain(self, domain):

        if isinstance(domain, Single):

            pos = domain.pid_particle_pair[1].position

            self.SinglesHistogram.bin(pos)
            self.NonMultisHistogram.bin(pos)

        elif isinstance(domain, Pair):

            pos1 = domain.pid_particle_pair1[1].position
            pos2 = domain.pid_particle_pair2[1].position

            for pos in [pos1, pos2]:
                self.PairsHistogram.bin(pos)
                self.NonMultisHistogram.bin(pos)

        elif isinstance(domain, Multi):

            for pid_particle_pair in domain.particle_container:

                pos = pid_particle_pair[1].position            
                self.MultisHistogram.bin(pos)

        else: 
            raise HistogramError( 'Cannot classify domain type for binning, domain = %s' % str(domain) )


    def normalize(self):
        
        for h in [self.SinglesHistogram, self.PairsHistogram, self.MultisHistogram, self.NonMultisHistogram]:

            try:

                h.normalize()
                h.check()

            except HistogramError as e:

                print "Warning: in %s: %s" % (self.__str__(), str(e))


    def __str__(self):

        return 'DomainsHistogramCollection \'' + str(self.name) + '\''



class HistogramError(Exception):
    pass

