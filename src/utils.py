import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.polynomial import polyfit
from eulerExact import riemann

def minmod( x , y , z ):
    prefactor = (1./4.) * np.abs( np.sign( x ) + np.sign ( y ) )
    minimum = np.minimum( np.abs( x ) , np.abs ( y ) )
    minimum = np.minimum( minimum , np.abs( z ) )
    return( prefactor * ( np.sign( x ) + np.sign( z ) ) * minimum )

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
