# Worker functions for fourierin

"""Worker function for fourierin.

Some day, these functions might be replaced by cython or other C code.

"""

import numpy as np

# 1d Fourier integral (FFT)
# -------------------------

def fourierin_1d_fft(f, a, b, c, d, r, s, midpoint = False):

    """Evaluate a 1d Fourier integral from function values at a regular spacing.

    The method implements the algorithm of Inverarity (2002), which is based on
    Bailey and Swarztrauber (1993), for the 1d case.

    Arguments
    ---------
    f : complex[:] or float[:]
        Function values
    a : float
        Lower integration boundary for t
    b : float
        Upper integration boundary for t
    c : float
        Lower limit of evaluation for w
    d : float
        Upper limit of evaluation for w
    r : float
        Factor for adjusting the constant factor of the FFT: 0 for $1/\sqrt{\pi}$,
        1 for no factor.
    s : float
        Factor for adjusting frequency: -1 or $- 2 \pi$ for the forward FFT, 1 or
        $2 \pi$ for the inverse FFT.
    midpoint : bool
        If True, function values are located on interval midpoints (default: False).

    Returns
    -------
    integral : complex[:]
        Fourier integral values for frequencies $c <= w < d$.

    Notes
    -----
    Notation/variable names are based on the notation in Inverarity (2002).

    References
    ----------
    Bailey, David H., and Paul N. Swarztrauber, 1994: A Fast Method for the Numerical Evaluation
       of Continuous Fourier and Laplace Transforms. SIAM Journal on Scientific Computing 15: 1105-1110.

    Inverarity, Gordon W., 2002: Fast Computation of Multidimensional Fourier Integrals. SIAM
       Journal on Scientific Computing, 24, pp. 645-651.

    """

    # Eqn. (4.2) for the handling of s

    c *= s  ;  d *= s

    # Eqns. (4.5) and (2.6)

    m = f.size  ;  beta = (b - a) / float(m)  ;  gamma = (d - c) / float(m)  ;  delta = 0.5*beta*gamma

    # Midpoints

    if midpoint:
        a += 0.5*beta  ;  b += 0.5*beta

    # Support arrays

    jdx = np.arange(0., float(m), 1.)
    t = a + beta * jdx  ;  w = c + gamma * jdx

    y     = np.zeros(2*m, dtype = complex)
    y[:m] = f * np.exp(1.j * jdx * (beta*c + delta*jdx))

    jdx   = np.arange(0., float(2*m), 1.)

    z     = np.empty(2*m, dtype = complex)
    z[:m] = np.exp(- 1.j * delta * jdx[:m]**2)
    z[m:] = np.exp(- 1.j * delta * (jdx[m:] - 2.*float(m))**2)

    # Convolution and inverse FFT - eq. (4.8)

    tmp = np.fft.ifft(np.fft.fft(y) * np.fft.fft(z))

    # Return result - rest eq. (4.8)

    return np.sqrt(np.abs(s) / (2.*np.pi)**(1.-r)) * beta * np.exp(1.j * (a*w + delta*jdx[:m]**2)) * tmp[:m]


# 1d complex Fourier integral, non-regular spacing
# ------------------------------------------------

def fourierin_1d_dft(f, a, b, w, r, s):

    """Evaluate a 1d Fourier integral using a DFT (for non-regular frequencies).

    The method implements a brute force evaluation of a Fourier integral for the 1d case.
    It may be useful if integral values are required for a small number of frequencies
    only.

    Arguments
    ---------
    f : complex[:]
        Function values
    a : float
        Lower integration boundary for t
    b : float
        Upper integration boundary for t
    w : float[:]
        Frequencies where to evaluate the Fourier integral
    r : float
        Factor for adjusting the constant factor of the FFT: 1 for no factor, -1 for 1/\sqrt(pi).
    s : float
        Factor for adjusting frequency: -1 or -2 pi for forward FFT, 1 or 2 pi for inverse FFT.

    Returns
    -------
    complex[:]
        Fourier integral values for frequencies w

    Notes
    -----
    Notation/variable names are based on the notation in Inverarity (2002).

    """

    m = f.size  ;  beta = (b - a) / float(m)  ;  k = w.size

    t   = np.linspace(a + beta/2., b - beta/2., num = m)
    fac = np.sqrt(np.abs(s) / (2.*np.pi)**(1-r)) * beta
    res = np.empty(k, dtype = complex)

    for i, wi in enumerate(w):
        res[i] = fac * np.sum(f * np.exp(1j * s * wi * t))

    return res