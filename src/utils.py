import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.polynomial import polyfit
from exact_riemann import riemann

def compute_l1_error( numerical , exact , deltaX ):
    diff = np.abs( numerical - exact )
    summ = diff.sum()
    return deltaX * summ

def fit_line( xdata , ydata ):
    # fit data to a line
    b, m = polyfit( xdata , ydata , 1 )
    return b , m

def line( x , b , m ):
    return m * x + b

def f(x,x0,sigma):
    return(np.where(abs(x-x0)<sigma , (1-((x-x0)/sigma)**2)**2 , 0))

def specific_internal_energy( rho , pressure , gamma=1.4 ):
    """ assumes ideal gas law """
    return( pressure / ( rho * ( gamma-1. ) ) )

def specific_enthalpy( rho , pressure , e ):
    return( 1 + e + (pressure / rho) )

def check_if_negative_pressures( pressures ):
    if isinstance(pressures,np.ndarray):
        if (pressures < 0.0).any():
            print("Warning, negative pressure encountered when computing sound speed")
    else:
        if pressures < 0.0:
            print("Negative pressure encountered when computing sound speed")
            print(pressures)

def lorentz_factor( v ):
    return( 1. / np.sqrt( 1. - v * v ) )
