# Main fourierin function

"""Main fourierin functions.

"""

import numpy as np

from .workers import fourierin_1d_fft, fourierin_1d_dft


# 1d Fourier integral
# -------------------

def fourierin_1d(f, lower_int, upper_int, lower_eval = None, upper_eval = None,
                 r = 1., s = -1., m = None, freqs = None,
                 midpoint = False, return_freqs = False, use_fft = True):

    """Evaluate a 1d Fourier integral with a complex integrand and regular spacing.

    The method implements Bailey & Swarztrauber (1993), based on Inverarity (2002),
    for the 1d case.

    Arguments
    ---------
    f : complex[:] or real[:]
        Function values
    lower_int : float
        Lower integration boundary for t
    upper_int : float
        Upper integration boundary for t
    lower_eval : float
        Lower limit of evaluation for w
    upper_eval : float
        Upper limit of evaluation for w
    r : float
        Factor for adjusting the constant factor of the FFT: 1 for no factor, -1 for 1/\sqrt(pi).
    s : float
        Factor for adjusting frequency: -1 or -2 pi for forward FFT, 1 or 2 pi for inverse FFT.
    m : int
        Number of values.
    freqs : float[:]
        Output frequencies; if given, a slow DFT algorithm will be used (default: None)
    midpoint : bool
        If True, function values are located on interval midpoints (default: False).
    return_freqs : bool
        If True, return frequencies (default: False).
    use_fft : bool
        If False, use a slow DFT rather then the FFT for the evaluation of the integral (default: True).

    Returns
    -------
    complex[:]
        Fourier integral values for frequencies c <= w < d

    Notes
    -----
    Notation/variable names are based on the notation in Inverarity (2002).

    """

    m = len(f)

    # FFT or brute-force DFT?

    if freqs is None:
        if lower_eval is None or upper_eval is None:
            raise ValueError("Unless freqs are provided, lower and upper evaluation bounds must be specified.")
    else:
        use_fft = False

    # Is f a function?

    if callable(f):
        if m is None:
            raise ValueError("If a function is given as input, the number of points to evaluate must be given.")
        beta = (upper_int - lower_int) / float(m)
        t    = np.linspace(lower_int + 0.5*beta, upper_int - 0.5*beta, num = m) # includes endpoint
        ft   = f(t)
        midpoint = True
    else:
        m  = len(f)
        ft = f

    # Evaluate the integral

    if use_fft:
        if return_freqs:
            gamma  = (upper_eval - lower_eval)/m
            freqs  = np.linspace(lower_eval, upper_eval - gamma, num = m)
            return fourierin_1d_fft(ft, lower_int, upper_int, lower_eval, upper_eval, r, s, midpoint = midpoint), freqs
        else:
            return fourierin_1d_fft(ft, lower_int, upper_int, lower_eval, upper_eval, r, s, midpoint = midpoint)
    else:
        return fourierin_1d_dft(ft, lower_int, upper_int, freqs, r, s)