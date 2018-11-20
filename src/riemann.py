
"""
solve the riemann problem
use a two shock approximation, and linearly interpolate
between the head and tail of the rarefaction
"""

URHO = 0
UMX = 1
UENER = 2

WRHO = 0
WV = 1
WP = 2

NVAR = 3


def riemann(wl, wr, gamma=1.4):
    """
    inputs:
    wl - input left primitives
    wr - input right primites
    """

    flux = np.zeros(NVAR)

    riemann_tolerance = 1.e-5
    nriemann = 15

    smallrho = 1.e-10
    smallp = 1.e-10
    smallv = 1.e-10

    rho_l = wl[WRHO]
    v_l = wl[WV]
    p_l = max(wl[WP],  smallp])

    rho_r = wr[WRHO]
    v_r = wr[WV]
    p_r= max(wr[WP],  smallp])

    # specific volume
    tau_l = 1. / max(rho_l, smallrho)
    tau_r = 1. / max(rho_r, smallrho)

    # sound speed
    c_l= np.sqrt(gamma * p_l * rho_l)
    c_r= np.sqrt(gamma * p_r * rho_r)

    # construct first guess for secant iteration by assuming
    # the nonlinear wave speed is equal to sound speed
    pstar1= p_r - p_l - c_r * (v_r - v_l)
    pstar1 = p_l + pstar1 * (c_l / (c_l + c_r))
    pstar1 = max(smallp, pstar1)

    # calculate nonlinear wave speeds for left and right
    # moving waves based on first guess for pressure
    # jump
    w_l1 = pstar1 + 0.5 * (gamma - 1.0) * (pstar1 + p_r)
    w_l1= np.sqrt(rho_l * abs(w_l1))

    w_r1= pstar1 + 0.5 * (gamma - 1.0) * (pstar1 + p_r)
    w_r1= np.sqrt(rho_r * abs(w_r1))

    # construct second guess for the pressure using
    # nonlinear wave speeds from the first guess
    pstar2 = p_r - p_l - w_r1 * (u_r - u_l)
    pstar2 = p_l + pstar2 * w_l1 / (w_l1 + w_r1)
    pstar2 = max(smallp, pstar2)

    # secant iteration
    # iterate until error falls below tolerance, or
    # until we hit nriem iterations

    has_converged = False

    for n in range(nriem):

        w_l = pstar2 + 0.5 * (gamma -1.0) * (pstar2 + p_l)
        w_l = np.sqrt( rho_l * abs( w_l ) )

        w_r = pstar2 + 0.5 * (gamma - 1.0) * (pstar2 + p_r )
        w_r = np.sqrt( rho_r * abs( w_r ) )

               ustar_l1 = u_l - (pstar1 - p_l)/w_l1
       ustar_r1 = u_r + (pstar1 - p_r)/w_r1
       ustar_l2 = u_l - (pstar2 - p_l)/w_l
       ustar_r2 = u_r + (pstar2 - p_r)/w_r

       delu1 = ustar_l1 - ustar_r1
       delu2 = ustar_l2 - ustar_r2

       scratch = delu2  - delu1

       if abs(pstar2 - pstar1) <= smallp:
           scratch = 0.

       if abs(scratch) < smallu:
          delu2 = 0.
          scratch = 1.

       # pressure at the "star" state -- using CG Eq. 18
       pstar = pstar2 - delu2*(pstar2 - pstar1)/scratch
       pstar = max(smallp, pstar)

       # check for convergence of iteration
       pres_err = abs(pstar - pstar2)/pstar
       if pres_err < riemann_tol:
          has_converged = True
          break

       # reset variables for next iteration
       pstar1 = pstar2
       pstar2 = pstar

       w_l1 = w_l
       w_r1 = w_r


    if not has_converged:
       print("Nonconvergence in subroutine rieman!")
       print("Pressure error = ", pres_err)
       print("pL = ", p_l, " pR = ", p_r)
       print("uL = ", u_l, " uR = ", u_r)
       print("cL = ", c_l, " c_r = ", c_r)
       print("Terminating execution")
       sys.exit("stopping")

    # end of secant iteration

    # calculate fluid velocity for the "star" state -- this comes from
    # the shock jump equations, Fry Eq. 68 and 69.  The ustar velocity
    # can be computed using either the jump eq. for a left moving or
    # right moving shock -- we use the average of the two.

    scratch = u_l - (pstar - p_l)/w_l
    scratch2 = u_r + (pstar - p_r)/w_r
    ustar = 0.5*(scratch + scratch2)

    if ustar < 0:
        ustar_sgn = -1.0
    elif ustar == 0.0:
        ustar_sgn = 0.0
    else:
        ustar_sgn = 1.0

    # decide which state is located at the zone iterface based on
    # the values of the wave speeds.  This is just saying that if
    # ustar > 0, then the state is U_L.  if ustar < 0, then the
    # state on the axis is U_R.
    scratch = 0.5*(1.0 + ustar_sgn)
    scratch2 = 0.5*(1.0 - ustar_sgn)

    ps = p_l*scratch + p_r*scratch2
    us = u_l*scratch + u_r*scratch2
    vs = tau_l*scratch + tau_r*scratch2

    rhos = 1.0/vs
    rhos = max(smlrho, rhos)

    vs = 1.0/rhos
    ws = w_l*scratch + w_r*scratch2
    ces = np.sqrt(gamma*ps*vs)

    # compute rhostar, using the shock jump condition (Fry Eq. 80)
    vstar = vs - (pstar - ps)/(ws*ws)
    rhostr = 1.0/ vstar
    cestar = np.sqrt(gamma*pstar*vstar)

    # compute some factors, Fry Eq. 81 and 82
    wes = ces - ustar_sgn*us
    westar = cestar - ustar_sgn*ustar

    scratch = ws*vs - ustar_sgn*us

    if pstar - ps >= 0.0:
       wes = scratch
       westar = scratch

    # compute correct state for rarefaction fan by linear interpolation
    scratch = max(wes - westar, wes + westar)
    scratch = max(scratch, smallu)

    scratch = (wes + westar)/scratch

    scratch = 0.5*(1.0 + scratch)
    scratch2 = 1.0 - scratch

    rhoav = scratch*rhostr + scratch2*rhos
    uav = scratch*ustar + scratch2*us
    pav = scratch*pstar + scratch2*ps

    if westar >= 0.0:
       rhoav  = rhostr
       uav = ustar
       pav = pstar

    if wes < 0.0:
       rhoav = rhos
       uav = us
       pav = ps

    # now compute the fluxes
    flux[URHO] = rhoav*uav
    flux[UMX] = rhoav*uav*uav + pav
    flux[UENER] = uav*(pav/(gamma - 1.0) + 0.5*rhoav*uav*uav + pav)

    return flux
