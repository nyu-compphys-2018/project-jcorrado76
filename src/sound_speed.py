from utils import specific_enthalpy, specific_internal_energy,\
check_if_negative_pressures
import numpy as np

def get_sound_speed(r , p, gamma ):
    """ 
    Inputs:
    r - (nparray) mass density 
    p - (nparray) pressure 
    gamma - (float) adiabatic constant 
    Returns:
    cs - (nparray) relativistic speed of sound
    """
    check_if_negative_pressures( p )
    e = specific_internal_energy( r , p , gamma=gamma )
    h = specific_enthalpy( r , p , e )
    return np.sqrt( ((gamma - 1.)/h) * (e + (p / r)) )
