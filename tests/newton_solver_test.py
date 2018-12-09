import numpy as np
from EulerSolver import specific_internal_energy,specific_enthalpy,fp,dfp,\
newton, lorentz_factor

r = 10.
v = 0.
p = 13.3
gamma = 1.4

print("Original primitives: \nr={}\nv={}\np={}".format(r,v,p))

lor = lorentz_factor( v )
e = specific_internal_energy( r , p , gamma )
h = specific_enthalpy( r , p , e )

print("Auxiliary variables: \nW={}\ne={}\nh={}".format( lor , e , h ) )

D = r * lor 
S = r * v * h * lor * lor # Sx
tau = r * h * lor**2 - p - lor * r  # tau

print( "Conserved variables: \nD={}\nS={}\ntau={}".format( D , S , tau ) )


p0 = 5
proot = newton( fp, dfp , p0 , D , S , tau , gamma)
vs = S / ( tau + proot + D )
lors = lorentz_factor( vs )
rs = D / lors

print("Recovered \nr={}\nv={}\np={}".format( np.allclose(r,rs),\
        np.allclose(v,vs),\
        np.allclose(p,proot)))

eps = 1e-3
print("Epsilon: " ,eps)
D = r * lor + eps
S = r * v * h * lor * lor + eps # Sx
tau = r * h * lor**2 - p - lor * r  + eps# tau

print("shifted D:",D)
print("shifted S:",S)
print("shifted tau:",tau)
print("shifted e:",e)
print("shifted h:",h)


p0 = 5
proot = newton( fp, dfp , p0 , D , S , tau , gamma)

vs = S / ( tau + proot + D )
lors = lorentz_factor( vs )
rs = D / lors

print("Recovered Shifted \nr={}\nv={}\np={}".format( np.allclose(r,rs),\
        np.allclose(v,vs),\
        np.allclose(p,proot)))

