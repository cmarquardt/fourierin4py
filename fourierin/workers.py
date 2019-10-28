# Worker functions for fourierin

"""Worker function for fourierin.

Some day, these functions might be replaced by cython or other C code.

"""

import numpy as np

# 1d complex Fourier integral
# ---------------------------

def fourierin_cmplx_1d(f, a, b, c, d, r, s):

    """Evaluate a 1d Fourier integral with a complex integrand and regular spacing.

    The method implements Bailey & Swarztrauber (1993), based on Inverarity (2002),
    for the 1d case.

    Arguments
    ---------
    f : complex[:]
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
        Factor for adjusting the constant factor of the FFT: 1 for no factor, -1 for 1/\sqrt(pi).
    s : float
        Factor for adjusting frequency: -1 or -2 pi for forward FFT, 1 or 2 pi for inverse FFT.

    Returns
    -------
    complex[:]
        Fourier integral values for frequencies c <= w < d

    Notes
    -----
    Notation/variable names are based on the notation in Inverarity (2002).

    """

    # Eqn. (4.2) for the handling of s

    c *= s  ;  d *= s

    # Eqns. (2.2) and (2.6)

    m = f.size  ;  beta = (b - a) / m  ;  gamma = (d - c) / m  ;  delta = 0.5*beta*gamma

    jdx = np.linspace(0., float(m-1), m)
    t = a + beta * jdx  ;  w = c + gamma * jdx

    y     = np.zeros(2*m)
    y[:m] = f * np.exp(1j * jdx * (beta*c + delta*jdx))

    jdx   = np.linspace(0., float(2*m - 1), 2*m)
    z     = np.empty(2*m)
    z[:m] = np.exp(1j * delta * jdx[:m]**2)
    z[m:] = np.exp(1j * delta * (jdx[m:] - 2*m)**2)

    # Convolution

    tmp = np.fft.ifft(np.fft.fft(y) * np.fft.fft(z))

    # Return result

    return np.sqrt(np.abs(s) / (2.*np.pi)**(1-r)) * beta * np.exp(1j * (a*w + delta*jdx[:m]**2)) * tmp[:m]


# 1d real Fourier integral
# ------------------------

def fourierin_real_1d(f, a, b, c, d, r, s):

    """Evaluate a 1d Fourier integral with a real integrand and regular spacing.

    The method implements Bailey & Swarztrauber (1993), based on Inverarity (2002),
    for the 1d case.

    Arguments
    ---------
    f : complex[:]
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
        Factor for adjusting the constant factor of the FFT: 1 for no factor, -1 for 1/\sqrt(pi).
    s : float
        Factor for adjusting frequency: -1 or -2 pi for forward FFT, 1 or 2 pi for inverse FFT.

    Returns
    -------
    complex[:]
        Fourier integral values for frequencies c <= w < d

    Notes
    -----
    Notation/variable names are based on the notation in Inverarity (2002).

    """

    raise NotImplementedError("This function is not yet implemented. Sorry!")


# 1d complex Fourier integral, non-regular spacing
# ------------------------------------------------

def fourierin_cmplx_1d_nonregular(f, a, b, c, d, r, s):

    """Evaluate a 1d Fourier integral with a complex integrand and non-regular spacing.

    The method implements Bailey & Swarztrauber (1993), based on Inverarity (2002),
    for the 1d case.

    Arguments
    ---------
    f : complex[:]
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
        Factor for adjusting the constant factor of the FFT: 1 for no factor, -1 for 1/\sqrt(pi).
    s : float
        Factor for adjusting frequency: -1 or -2 pi for forward FFT, 1 or 2 pi for inverse FFT.

    Returns
    -------
    complex[:]
        Fourier integral values for frequencies c <= w < d

    Notes
    -----
    Notation/variable names are based on the notation in Inverarity (2002).

    """

    raise NotImplementedError("This function is not yet implemented. Sorry!")


# 1d real Fourier integral, non-regular spacing
# ---------------------------------------------

def fourierin_real_1d_nonregular(f, a, b, c, d, r, s):

    """Evaluate a 1d Fourier integral with a real integrand and non-regular spacing.

    The method implements Bailey & Swarztrauber (1993), based on Inverarity (2002),
    for the 1d case.

    Arguments
    ---------
    f : complex[:]
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
        Factor for adjusting the constant factor of the FFT: 1 for no factor, -1 for 1/\sqrt(pi).
    s : float
        Factor for adjusting frequency: -1 or -2 pi for forward FFT, 1 or 2 pi for inverse FFT.

    Returns
    -------
    complex[:]
        Fourier integral values for frequencies c <= w < d

    Notes
    -----
    Notation/variable names are based on the notation in Inverarity (2002).

    """

    raise NotImplementedError("This function is not yet implemented. Sorry!")