# Main fourierin function

"""Main fourierin functions.

"""

import numpy as np

from .workers import fourierin_1d_fft, fourierin_1d_dft


# 1d Fourier integral
# -------------------

def fourierin_1d(f, lower_int, upper_int, lower_eval = None, upper_eval = None,
                 const_adj = 1., freq_adj = -1., resolution = None,
                 omegas = None, use_fft = True):

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
    const_adj : float
        Factor for adjusting the constant factor of the FFT: 1 for no factor, -1 for 1/\sqrt(pi).
    freq_adj : float
        Factor for adjusting frequency: -1 or -2 pi for forward FFT, 1 or 2 pi for inverse FFT.
    resolution : int
        Number of values.

    Returns
    -------
    complex[:] or real[:]
        Fourier integral values for frequencies c <= w < d

    Notes
    -----
    Notation/variable names are based on the notation in Inverarity (2002).

    """

    # TODO: support f as a function as in the R code

    resolution = len(f)

    # Frequencies and their increment

    if omegas is None:
        if lower_eval is None or upper_eval is None:
            raise ValueError("Unless omegas are provided, lower and upper evaluation bounds must be specified.")

        gamma  = (upper_eval - lower_eval)/resolution
        omegas = np.linspace(lower_eval, upper_eval - gamma, num = resolution)
    else:
        use_fft = False

    # Function evaluation (not yet implemented)

    ### If f is the function, it needs to be evaluated in
    ### the time domain values.
    #if (is.function(f)) {

    #    del <- (b - a)/resol # Increment in the time domain.
    #    t <- seq(a + del/2, b - del/2,
    #             length.out = resol)    # Freq. dom. vector.
    #    f_t <- f(t)                     # Function values
    #    ## Rutinary check
    #    if(is.null(f_t)) stop("Function f is null.")
    #} else {
    #    f_t <- f
    #}

    f_t = f

    # Evaluate the integral

    if not use_fft:
        return fourierin_1d_nonregular(f_t, lower_int, upper_int, lower_eval, upper_eval, const_adj, freq_adj)
    else:
        return fourierin_1d(f_t, lower_int, upper_int, lower_eval, upper_eval, const_adj, freq_adj)