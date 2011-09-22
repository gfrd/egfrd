#!/usr/bin/env python

__author__    = 'Koichi Takahashi <shafi@e-cell.org>'
__license__   = 'GPL'
__copyright__ = 'Copyright The Molecular Sciences Institute 2006-2007'


import unittest

import _greens_functions as mod

import math
import numpy


class FreeFunctionsTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_int_p_theta_free_is_ip_theta_free(self):

        import scipy.integrate

        D = 1e-12
        t = 1e-5
        sigma = 1e-9
        r0 = 1e-9
        r = r0
        kf = 1e-18
        
        ip = mod.ip_theta_free(0.0, r, r0, t, D)
        self.assertEqual(0.0, ip)
        
        resolution = 10
        for i in range(1, resolution):
            theta = i * numpy.pi / resolution 
            ip = mod.ip_theta_free(theta, r, r0, t, D)
            result = scipy.integrate.quad(mod.p_theta_free, 0.0, theta,
                                          args=(r, r0, t, D))
            np = result[0]
            self.assertAlmostEqual(0.0, (np-ip)/ip)


    def test_int_p_irr_is_p_survival_irr(self):

        import scipy.integrate

        D = 1e-12
        t = 1e-5
        sigma = 1e-9
        r0 = 1e-9
        kf = 1e-18
        

        for i in range(1, 20):
            S = mod.p_survival_irr(t, r0 * i, kf, D, sigma)
            result = scipy.integrate.quad(mod.p_irr, sigma, sigma * 1e3,
                                          args=(t, r0 * i, kf, D, sigma))
            ip = result[0]
            self.failIf(ip == 0)
            self.assertAlmostEqual(0.0, (S-ip)/ip)


    def test_int_g_bd_is_I_bd(self):

        import scipy.integrate
        import math

        D = 1e-12
        t = 1e-11
        sigma = 2e-9
        r0 = 2e-9
        v = 0.01

        ibd = mod.I_bd_3D(sigma, t, D)
        ibd2 = mod.I_bd_1D(sigma, t, D, v)
        #print ibd, ibd2

        result = scipy.integrate.quad(mod.g_bd_3D, sigma, 
                                      sigma + 6 * math.sqrt(6 * D * t),
                                      args=(sigma, t, D))

        result2 = scipy.integrate.quad(mod.g_bd_1D, sigma, 
                                      sigma + 100 * math.sqrt(2 * D * t),
                                      args=(sigma, t, D, v))
 
        igbd = result[0]
        igbd2 = result2[0]

        #print igbd2, igbd2
        self.failIf(ibd == 0)
        self.failIf(ibd2 == 0)
        self.assertAlmostEqual(0.0, (ibd-igbd)/ibd)
        #self.assertAlmostEqual(0.0, (ibd2-igbd2)/ibd2)


    def test_int_g_bd_is_I_bd_smallt(self):

        import scipy.integrate

        D = 1e-12
        t = 1e-20
        sigma = 1e-8
        r0 = 1e-9
        v = 0.01

        ibd = mod.I_bd_3D(sigma, t, D)
        ibd2 = mod.I_bd_1D(sigma, t, D, v)
        #print ibd, ibd2

        result = scipy.integrate.quad(mod.g_bd_3D, sigma, sigma + 
                                      6 * math.sqrt(6 * D * t),
                                      args=(sigma, t, D))

        result2 = scipy.integrate.quad(mod.g_bd_1D, sigma, 
                                      sigma + 10 * math.sqrt(2 * D * t),
                                      args=(sigma, t, D, v))

        igbd = result[0]
        igbd2 = result2[0]

        #print igbd, igbd2
        self.failIf(ibd == 0)
        self.failIf(ibd2 == 0)
        self.assertAlmostEqual(0.0, (ibd-igbd)/ibd)
        self.assertAlmostEqual(0.0, (ibd2-igbd2)/ibd2)


    def test_I_bd_r_large_is_I_bd(self):

        D = 1e-12
        t = 1e-11
        sigma = 2e-9
        r0 = 5e-9
        v = 0.01

        ibd = mod.I_bd_3D(sigma, t, D)
        ibdr = mod.I_bd_r_3D(sigma + 6 * math.sqrt(6 * D * t), sigma, t, D)

        ibd2 = mod.I_bd_1D(sigma, t, D, v)
        ibdr2 = mod.I_bd_r_1D(sigma + 10 * math.sqrt(2 * D * t), sigma, t, D, v)
        #print ibd, ibdr
        #print ibd2, ibdr2

        self.assertAlmostEqual(0.0, (ibd-ibdr)/ibd)
        self.assertAlmostEqual(0.0, (ibd2-ibdr2)/ibd2)


    def test_int_g_bd_is_I_bd_r(self):

        import scipy.integrate
        import math

        D = 1e-12
        t = 1e-11
        sigma = 2e-9
        v = 0.01

        r_max = 6 * math.sqrt(6 * D * t)

        ibd = mod.I_bd_r_3D(sigma, sigma, t, D)
        self.failIf(ibd != 0.0)
        ibd = mod.I_bd_r_1D(sigma, sigma, t, D, v)
        self.failIf(ibd != 0.0)

        N = 20
        for i in range(1, N):
            r = sigma + r_max / N * i
            ibd = mod.I_bd_r_3D(r, sigma, t, D)
            ibd2 = mod.I_bd_r_1D(r, sigma, t, D, v)

            result = scipy.integrate.quad(mod.g_bd_3D, sigma, r,
                                          args=(sigma, t, D))

            result2 = scipy.integrate.quad(mod.g_bd_1D, sigma, r,
                                          args=(sigma, t, D, v))

            igbd = result[0]
            igbd2 = result2[0]

            self.failIf(ibd == 0)
            self.assertAlmostEqual(0.0, (ibd-igbd)/ibd)
            self.failIf(ibd2 == 0)
            self.assertAlmostEqual(0.0, (ibd2-igbd2)/ibd2)

    def test_drawR_gbd(self):

        import scipy.integrate
        import math

        D = 1e-12
        t = 1e-11
        sigma = 2e-9
        v = 0.01

        r = mod.drawR_gbd_3D(0.0, sigma, t, D)
        self.assertEqual(r, sigma)

        r = mod.drawR_gbd_3D(0.5, sigma, t, D)
        self.failIf(r <= sigma)
        #print 'rr', r

        r = mod.drawR_gbd_3D(1.0, sigma, t, D)
        self.failIf(r <= sigma)

        r = mod.drawR_gbd_1D(0.0, sigma, t, D, v)
        #print 'rr@0.0', r
        self.assertEqual(r, sigma)

        r = mod.drawR_gbd_1D(0.5, sigma, t, D, v)
        #print 'rr@0.5', r
        self.failIf(r <= sigma)

        r = mod.drawR_gbd_1D(1.0, sigma, t, D, v)
        #print 'rr@1.0', r
        self.failIf(r <= sigma)


    def test_p_reaction_irr_t_inf(self):
        
        D = 1e-12
        t = numpy.inf
        sigma = 1e-8
        r0 = 1.1e-8
        kf = 1e-16
        kr = 10
        kD = 4 * numpy.pi * sigma * D

        alpha = (1 + (kr / kD)) * math.sqrt(D) / sigma

        pr = mod.p_reaction_irr(t, r0, kf, D, sigma, alpha, kD)
        prinf = mod.p_reaction_irr_t_inf(r0, kf, sigma, kD)

        #print pr, prinf

        self.assertAlmostEqual(pr, prinf)


        
if __name__ == "__main__":
    unittest.main()
