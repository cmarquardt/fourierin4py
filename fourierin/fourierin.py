# Main fourierin function

"""Main fourierin functions.

"""

import numpy as np

from .workers import fourierin_cmplx_1d, fourierin_cmplx_1d_nonregular, fourierin_real_1d, fourierin_real_1d_nonregular


"""
fourierin_1d <- function(f, lower_int, upper_int,
                         lower_eval = NULL, upper_eval = NULL,
                         const_adj, freq_adj, resolution = NULL,
                         eval_grid = NULL, use_fft = TRUE) {
    ## Condensed notation
    a <- lower_int
    b <- upper_int
    c <- lower_eval
    d <- upper_eval
    r <- const_adj
    s <- freq_adj
    resol <- resolution
    w <- eval_grid

## fourierin_1d <- function(f, a, b, c = NULL, d = NULL,
##                          r, s, resol = NULL, w = NULL,
##                          use_fft = TRUE) {
    ## Flag to determine wheter a list or a vector will be returned
    w_given <- !is.null(w)

    ## If function values are provided, then the resolution
    ## is the length of the vector of values.
    if (!is.function(f)) resol <- length(f)

    ## Increment in the frequency domain.
    gam <- (d - c)/resol

    ## Freq. dom. vector. If w is provided, FFT will NOT be used.
    if (is.null(w)) {
        if (is.null(c) | is.null(d)) {
            stop("c and d must be provided.")
        }
        w <- seq(c, d - gam, length.out = resol)
    } else {
        use_fft <- FALSE
    }

    ## If f is the function, it needs to be evaluated in
    ## the time domain values.
    if (is.function(f)) {

        del <- (b - a)/resol # Increment in the time domain.
        t <- seq(a + del/2, b - del/2,
                 length.out = resol)    # Freq. dom. vector.
        f_t <- f(t)                     # Function values
        ## Rutinary check
        if(is.null(f_t)) stop("Function f is null.")
    } else {
        f_t <- f
    }

    if (!use_fft) {
        out <- switch(is.complex(f_t) + 1,
                      fourierin_1d_nonregular_cpp(f_t, a, b, w,
                                                  resol, r, s),
                      fourierin_cx_1d_nonregular_cpp(f_t,
					 a, b, w, resol, r, s))
    } else {
        out <- switch(is.complex(f_t) + 1,
                      fourierin_1d_cpp(f_t, a, b, c, d, r, s),
                      fourierin_cx_1d_cpp(f_t, a, b, c, d, r, s))
    }

    ## If w is given, return only the values of the integral,
    ## otherwise alse return w.
    if(w_given) return(out)

    return(list(w = w,                  # Return list.
                values = drop(out)))
}"""

# 1d complex Fourier integral
# ---------------------------

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
        if f.dtype == np.dtype(complex):
            return fourierin_cmplx_1d_nonregular(f_t, lower_int, upper_int, lower_eval, upper_eval, const_adj, freq_adj)
        else:
            return fourierin_real_1d_nonregular(f_t, lower_int, upper_int, lower_eval, upper_eval, const_adj, freq_adj)
    else:
        if f.dtype == np.dtype(complex):
            return fourierin_cmplx_1d(f_t, lower_int, upper_int, lower_eval, upper_eval, const_adj, freq_adj)
        else:
            return fourierin_cmplx_1d(f_t, lower_int, upper_int, lower_eval, upper_eval, const_adj, freq_adj)