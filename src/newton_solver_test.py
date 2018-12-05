import numpy as np
from EulerSolver import specific_internal_energy,specific_enthalpy,fp,dfp,\
newton, lorentz_factor

r = 10.
v = 0.
p = 13.3
gamma = 1.4

lor = lorentz_factor( v )
e = specific_internal_energy( r , p , gamma )
h = specific_enthalpy( r , p , e )


eps = 1e-3


D = r * lor + eps
S = r * v * h * lor * lor + eps # Sx
tau = r * h * lor**2 - p - lor * r  + eps# tau

print("D:",D)
print("S:",S)
print("tau:",tau)
print("e:",e)
print("h:",h)


p0 = 5
proot = newton( fp, dfp , p0 , D , S , tau , gamma)

vs = S / ( tau + proot + D )
lors = lorentz_factor( vs )
rs = D / lors

print("r:",rs)
print("v:",vs)
print("p:",proot)
print( np.allclose( r , rs ) )
print( np.allclose( v , vs ) )
print( np.allclose( p , proot ) )
