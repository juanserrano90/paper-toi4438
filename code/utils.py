import numpy as np
import celerite2

# this functions were made to compute the LOO log predictive probability for a Gaussian process regression model
# for the case 

# manually compute the LOO log probability density

def stabilize_matrix(a, jitter=1e-6):
    """
    Add a small jitter to the diagonal of a matrix to stabilize it.
    """
    a[np.arange(a.shape[0]), np.arange(a.shape[1])] += jitter
    # a = a + jitter
    return a
def add_diagonal(a, diag):
    """
    Add a diagonal to a matrix.
    """
    a[np.arange(a.shape[0]), np.arange(a.shape[1])] += diag
    return a

def lag_matrix(y, x):  # y is the one with one missing value
    """
    Compute the lag matrix: the difference between `y` and `x`.
    """
    matrix = y[:, None] - x[None, :]
    return matrix

# def sho_term(sigma, rho, Q=1.0 / 3):
#     """Create a SHOTerm kernel with given parameters."""
#     return celerite2.terms.SHOTerm(sigma=sigma, rho=rho, Q=Q)

def sho_term(rho, tau, sigma):
    """Create a SHOTerm kernel with given parameters."""
    return celerite2.terms.SHOTerm(rho=rho, tau=tau, sigma=sigma)

def rotation_term(sigma, period, Q0, dQ, f):
    """Create a Rotation term kernel with given parameters."""
    return celerite2.terms.RotationTerm(sigma=sigma, period=period, Q0=Q0, dQ=dQ, f=f)

def logp(indice, x, y, diag, sigma, term=None, tau=None, rho=None, period=None, Q0=None, dQ=None, f=None, verbose=False):
    """Compute the LOO log predictive probability for a single point."""
    # instanciate the kernel
    if term == 'sho':
        # kernel = sho_term(sigma, rho)
        kernel = sho_term(rho, tau, sigma)
    elif term == 'rotation':
        kernel = rotation_term(sigma, period, Q0, dQ, f)
    else:
        raise ValueError("term must be 'sho' or 'rotation'")

    # compute the intermediate matrices
    # m1 = kernel.get_value(lag_matrix(np.array([x[indice]]), np.delete(x,indice)))  # m1 has small values, makes sense for pairs of points too separated, ~1e-6
    # m2 = kernel.get_value(lag_matrix(np.delete(x,indice), np.delete(x,indice)))
    # m2 = add_diagonal(m2, np.delete(diag, indice))
    # m3 = kernel.get_value(lag_matrix(np.array([x[indice]]), np.array([x[indice]])))    # the input to get_value here is always 0, because there's only one test point in  K(x*, x*)
    # m4 = kernel.get_value(lag_matrix(np.delete(x,indice), np.array([x[indice]])))

    # compute the mean and covariance of the predictive distribution for the left-out point
    # mean = np.dot(m1, np.linalg.solve(m2, np.delete(y, indice)))
    # cov = m3 - np.dot(m1, np.linalg.solve(m2, m4))

    # alternative way
    m5 = kernel.get_value(lag_matrix(x, x))
    m5 = add_diagonal(m5, diag)

    # mean = y[indice] - np.linalg.solve(m5, y)[indice] / np.linalg.solve(m5, np.identity(len(x)))[indice, indice]
    # cov = 1 / np.linalg.solve(m5, np.identity(len(x)))[indice, indice]

    m5_inv = np.linalg.solve(m5, np.identity(len(x)))
    mean = y[indice] - np.dot(m5_inv, y)[indice] / m5_inv[indice, indice]
    cov = 1 / m5_inv[indice, indice]

    # compute the log predictive probability
    term1 = -0.5*np.log(cov)
    term2 = -0.5*(y[indice]-mean)**2/cov 
    term3 = -0.5*np.log(2*np.pi)
    
    if verbose:
        print('indice:', indice, 'cov:', cov, 'mean:', mean, "term1:", term1, "term2:", term2, "term3:", term3, "diag:", diag[indice])
    return term1 + term2 + term3

def loo_cv(x, y, diag, sigma, term=None, rho=None, tau=None, period=None, Q0=None, dQ=None, f=None, verbose=False):
    """Compute the full LOO-CV log predictive probability.
    term (str): 'sho' or 'rotation'
    """
    logp_values = [logp(i, x, y, diag, sigma, term=term, rho=rho, tau=tau, period=period, Q0=Q0, dQ=dQ, f=f, verbose=verbose) for i in range(len(x))]

    loo = np.sum(logp_values)
    loo_se = (len(x)*np.var(logp_values))**0.5
    return loo, loo_se

def log_likelihood(x, y, diag, term, sigma, rho=None, period=None, Q0=None, dQ=None, f=None):
    """log-Likelihood function for N data points yn at points tn
    with variance diag and kernel term.
    this does the same as the compute_log_likelihood function in pymc"""
    if term == 'sho':
        kernel = sho_term(sigma, rho)
    elif term == 'rotation':
        kernel = rotation_term(sigma, period, Q0, dQ, f)
    else:
        raise ValueError("term must be 'sho' or 'rotation'")
    
    m1 = kernel.get_value(lag_matrix(x, x))
    m1 = add_diagonal(m1, diag)

    term1 = -0.5*np.dot(y.T, np.linalg.solve(m1, y))
    term2 = -0.5*np.log(np.linalg.det(m1))
    term3 = -0.5*len(x)*np.log(2*np.pi)

    return term1 + term2 + term3

def tsm(teq, r_p, m_p, m_j, r_star,):
    ''' Function to calculate the Transmission Spectroscopy Metric (TSM) for exoplanets 
        as defined by Kempton et al. (2018) 
        TSM = (Scale factor) x (r_p**3 Teq)/(m_p r_star**2) x 10^(-m_J/5)

        where:
        teq: Equilibrium temperature of the planet in Kelvin (calculated assuming zero albedo and full heat redistribution)
        r_p: Radius of the planet in Earth radii
        m_p: Mass of the planet in Earth masses
        m_J: J-band magnitude of the host star, chosen as a filter that is close to the middle of the NIRISS bandpass
    '''
    if r_p < 1.5:
        scale_factor = 0.190
    elif 1.5 <= r_p < 2.75:
        scale_factor = 1.26
    elif 2.75 <= r_p < 4.0:
        scale_factor = 1.28
    elif 4.0 <= r_p < 10.0:
        scale_factor = 1.15
    else:
        raise ValueError("Planet radius out of range for TSM calculation")
    
    tsm_value = scale_factor * (r_p**3 * teq) / (m_p * r_star**2) * 10**(-m_j / 5)

    return tsm_value

def esm(teq, r_p, t_star, r_star, m_k):
    ''' Function to calculate the Emission Spectroscopy Metric (ESM) for exoplanets 
        as defined by Kempton et al. (2018) 
        ESM = 4.29 x 10^6 x (B_7.5(Tday)/B_7.5(Tstar)) x (r_p/r_star)^2 x 10^(-m_K/5)

        where:
        Tday: the planet's dayside temperature in Kelvin, we calculate it as 1.10 x Teq
        B_7.5: Planck function evaluated at 7.5 microns for a given temperature
        r_p: Radius of the planet in Earth radii
        r_star: Radius of the host star in Solar radii
        m_K: the apparent mag of the star in the K band
    '''
    # Planck function at 7.5 microns
    def planck_7_5(temp):
        h = 6.62607015e-34  # Planck constant in J*s
        c = 2.99792458e8     # Speed of light in m/s
        k = 1.380649e-23  # Boltzmann constant in J/K
        wavelength = 7.5e-6  # Wavelength in meters

        B = (2*h*c**2) / (wavelength**5) * 1 / (np.exp((h*c) / (wavelength*k*temp)) - 1)
        return B
    
    T_day = 1.10 * teq
    B_day = planck_7_5(T_day)
    B_tstar = planck_7_5(t_star)

    esm_value = 4.29e6 * (B_day / B_tstar) * (r_p / r_star)**2 * 10**(-m_k / 5)
    
    return esm_value