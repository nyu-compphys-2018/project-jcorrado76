from utils import specific_internal_energy, specific_enthalpy, lorentz_factor
from numpy import zeros
def prim_to_cons( W ,gamma):
    """ compute relativistic conserved variables """
    U = zeros(W.shape)
    r = W[0,:]
    v = W[1,:]
    p = W[2,:]
    e = specific_internal_energy( r , p , gamma )
    h = specific_enthalpy( r , p , e )
    lorentz = lorentz_factor(v)
    U[0,:] = lorentz * r  # D
    U[1,:] = r * v * h * lorentz * lorentz # Sx
    U[2,:] = r * h * lorentz * lorentz - p - lorentz * r # tau
    return U
