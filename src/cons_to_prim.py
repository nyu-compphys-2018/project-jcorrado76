from newton import newton, fp, dfp
from numpy import zeros, where , fabs
from utils import lorentz_factor

def cons_to_prim( U, gamma ):
    """ perform a recovery of the primitive variables """
    W = zeros(U.shape)
    D = U[0,:]
    S = U[1,:]
    tau = U[2,:]
    ps = zeros(U.shape[1])
    p0 = 5
    for i in range(U.shape[1]):
        proot = newton( fp, dfp , p0 , D[i] , S[i],tau[i], gamma)
        ps[i] = proot

    vs = S / ( tau + ps + D )
    lors = lorentz_factor( vs )
    rs = D / lors

    if (ps<0).any():
        print("Negative pressure")
        print(ps)
        ps = where(ps<0,fabs(S-tau-D),ps)
    if ( vs >= 1).any():
        print( "v greater than c")
        print(vs[vs>1])
        vs = where( vs>1 , 1e-6 , vs )

    W[2,:] = ps[:]
    W[1,:] = vs[:]
    W[0,:] = rs[:]

    return W
