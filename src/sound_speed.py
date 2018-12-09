from utils import specific_enthalpy, specific_internal_energy,\
check_if_negative_pressures
import numpy as np

def get_sound_speed(r , p, gamma ):
    """ get relativistic sound speed """
    check_if_negative_pressures( p )
    e = specific_internal_energy( r , p , gamma=gamma )
    h = specific_enthalpy( r , p , e )
    return np.sqrt( ((gamma - 1.)/h) * (e + (p / r)) )
