#------------------------------------------------------------------------------
# 1. Imports
#------------------------------------------------------------------------------

import numpy     as np
import fourierin as fi


#------------------------------------------------------------------------------
# 2. Useful things
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# 2. 1d tests - fft
#------------------------------------------------------------------------------

class Test1dGaussian(object):

    def setup_method(self):

        n     = 256
        self.t_min = -10.  ;  self.t_max = 10.  ;  delta_t = (self.t_max - self.t_min) / float(n)
        self.w_min =  -5.  ;  self.w_max =  5.  ;  delta_w = (self.w_min - self.w_max) / float(n)

        self.t = np.linspace(self.t_min, self.t_max - delta_t, num = n)  # !!!NOTE THE CALCULATION OF TIME STAMPS
        self.w = np.linspace(self.w_min, self.w_max - delta_w, num = n)  #    AND FREQUENCIES!!!

        self.f = np.exp(- self.t ** 2 / 2.) / np.sqrt(2. * np.pi)
        self.g = np.exp(- self.w ** 2 / 2.) / np.sqrt(2. * np.pi)


    def test_gaussian_1d_fft(self):

        g = fi.fourierin_1d(self.f, self.t_min, self.t_max, self.w_min, self.w_max, r = 0., s = -1.)

        np.allclose(g, self.g, rtol = 1.e-10, atol = 1.e-10)


    def test_gaussian_1d_fft_freqs(self):

        g, w = fi.fourierin_1d(self.f, self.t_min, self.t_max, self.w_min, self.w_max, r = 0., s = -1., return_freqs = True)

        np.allclose(g, self.g, rtol = 1.e-10, atol = 1.e-10)
        np.allclose(w, self.w, rtol = 1.e-10, atol = 1.e-10)


    def test_gaussian_1d_dft(self):

        idx = np.array([12, 87, 116, 214, 252])
        w   = self.w[idx]

        g = fi.fourierin_1d(self.f, self.t_min, self.t_max, freqs = w, r = 0., s = -1.)

        np.allclose(g, self.g[idx], rtol = 1.e-10, atol = 1.e-10)