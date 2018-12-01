import numpy as np
def minmod( a , b ):
    if abs(a) < abs(b) and a * b > 0.0:
        return a
    elif abs(b) < abs(a) and a * b > 0.0:
        return b
    else:
        return 0.0
def maxmod( a , b ):
    if abs(a) > abs(b) and a * b > 0.0:
        return a
    elif abs(b) < abs(a) and a * b > 0.0:
        return b
    else:
        return 0.0
def shifted(s, x):
    # returns slice shifted by x
    return slice(s.start + x, s.stop + x)


class State_Reconstructor( object ):
    def __init__(self, grid,type='minmod'):
        self.grid=grid
        slope = g.get_scratch_array()

    def compute_states(self,dt):
        reconstruction = slice(g.ilo-1,g.ihi+2)
        if self.type == "godunov":
            # piecewise constant slope
            slope[:,:] = 0.0
        elif self.type == "centered":
            # central difference
            slope[:,reconstruction] = 0.5 * (g.U[:,shifted(reconstruction,1)] - g.U[:,shifted(reconstruction,-1)] ) / g.dx
        elif self.type == "minmod":
            # minmod
            slope[:,reconstruction] = minmod( (g.U[:,reconstruction] - g.U[:,shifted(reconstruction,-1)])/g.dx,
                               (g.U[:,shifted(reconstruction,1)] - g.U[reconstruction])/g.dx )
        elif self.type == "MC":
            # MC limiter
            slope[reconstruction] = minmod(minmod( 2.0*(g.U[:,reconstruction] - g.U[:,shifted(reconstruction,-1)])/g.dx,
                                      2.0*(g.U[:,shifted(reconstruction,1)] - g.U[:,reconstruction])/g.dx ),
                              0.5*(g.U[:,shifted(reconstruction,1)] - g.U[:,shifted(reconstruction,-1)])/g.dx)
        elif self.type == "superbee":
            # superbee limiter
            A = minmod( (g.U[:,shifted(reconstruction,1)] - g.U[:,reconstruction])/g.dx,
                        2.0*(g.U[:,reconstruction] - g.U[:,shifted(reconstruction,-1)])/g.dx )

            B = minmod( (g.U[:,reconstruction] - g.U[:,shifted(reconstruction,-1)])/g.dx,
                        2.0*(g.U[shifted(reconstruction,1)] - g.U[:,reconstruction])/g.dx )

            slope[i] = maxmod(A, B)

        UL = g.get_scratch_array()
        UR = g.get_scratch_array()

        UR[ : , g.ilo-1 : g.ihi + 3 ] = U[ : , g.ilo-1 : g.ihi + 3 ] - \
                0.5 * ( 1.0 + U[ : , g.ilo-1 : g.ihi+3 ] * dt / self.grid.dx ) * \
                ldeltau[ : , g.ilo-1 : g.ihi+3 ]

        UL[ : , g.ilo : g.ihi+3 ] = U[ : , g.ilo-1 : g.ihi+2 ] +\
                0.5 * ( 1.0 - U[ : , g.ilo-1 : g.ihi+2 ] * dt / self.grid.dx ) * \
                ldeltau[ : , g.ilo-1 : g.ihi+2 ]

        return UL, UR

        def MC( self, dt ):
            """ reconstruct left and right interface states """
            g = self.grid

            # compute piecewise linear slopes -- 2nd order MC limiter
            # pick a range of cells that includes 1 ghost cell on either side
            ie = g.ihi + 1

            U = g.U

            dc = g.get_scratch_array()
            dl = g.get_scratch_array()
            dr = g.get_scratch_array()

            dc[ : , g.ilo-1:ie+1] = 0.5* ( U[ : , g.ilo : ie+2 ] - U[ : , g.ilo-2 : ie ] )
            dl[ : , g.ilo-1:ie+1] = U[ : , g.ilo-1+1 : ie+2 ] - U[ : , g.ilo-1 : ie+1 ]
            dr[ : , g.ilo-1:ie+1] = U[ : , g.ilo-1 : ie+1 ] - U[ : , g.ilo-2 : ie ]

            # minmod
            d1 = 2.0 * np.where( np.fabs( dl ) < np.fabs( dr ) , dl , dr )
            d2 = np.where( np.fabs( dc ) < np.fabs( d1 ) , dc , d1 )
            ldeltau = np.where( dl * dr > 0.0 , d2 , 0.0 )

            # now interface states. there is one more interface than zones
            ul = g.get_scratch_array()
            ur = g.get_scratch_array()

            ur[ : , g.ilo-1 : ie + 2 ] = U[ : , g.ilo-1 : ie + 2 ] - \
                    0.5 * ( 1.0 + U[ : , g.ilo-1 : ie+2 ] * dt / self.grid.dx ) * \
                    ldeltau[ : , g.ilo-1 : ie+2 ]

            ul[ : , g.ilo : ie+2 ] = U[ : , g.ilo-1 : ie+1 ] +\
                    0.5 * ( 1.0 - U[ : , g.ilo-1 : ie+1 ] * dt / self.grid.dx ) * \
                    ldeltau[ : , g.ilo-1 : ie+1 ]

            return ul , ur

        def textbook_minmod( self ):
            """ minmod slope limiter """
            g = self.grid
            UL = g.get_scratch_array()
            UR = g.get_scratch_array()
            Ujm1 = U[:,shifted(g.physical,-1)]
            Uj   = U[:,g.physical]
            Ujp1 = U[:,shifted(g.physical,1)]
            sp = Ujp1 - Uj
            sm = Uj - Ujm1
            ssp = np.sign(sp)
            ssm = np.sign(sm)
            asp = np.abs( sp )
            asm = np.abs( sm )
            dU = 0.25 * ( ssp + ssm ) * np.minimum( asp , asm )
            Ujp = Uj + dU
            Ujm = Uj - dU
            UL[:,g.physical] = Ujp
            UR[:,shifted(g.physical,-1)] = Ujm
            return UL, UR
