from numpy import sqrt
from utils import lorentz_factor,specific_internal_energy,specific_enthalpy
def fp( p , D , S , tau , gamma=1.4):
    vstar = S / ( tau + p + D )
    wstar = 1. / sqrt( 1 - vstar * vstar )
    rstar = D / wstar
    estar = ( tau + D * ( 1. - wstar ) + ( 1 - wstar * wstar ) * p ) / ( D * wstar )
    return (( gamma - 1.) * rstar * estar - p)

def dfp(p , D , S , tau , gamma=1.4):
    """ approximate derivative of the above expression """
    vstar = S / ( tau + p + D )
    lor = lorentz_factor( vstar )
    r = D / lor
    e = specific_internal_energy( r , p , gamma=gamma )
    h = specific_enthalpy( r , p , e )
    cs2 = ((gamma - 1.)/h) * (e + (p / r))
    return(vstar * vstar * cs2 - 1)

def newton( func , fprime , x0 , *args, tol=1e-12 ):
    delta = abs( 0.0 - func(x0 , *args) )
    while delta > tol:
        x0 = x0 - func( x0 , *args) / fprime( x0 , *args)
        delta = abs( 0.0 - func(x0, *args) )
    root = x0
    return root
