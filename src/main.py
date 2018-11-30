def main():
    # PARAMETERS
    xmin = 0.0
    xmax = 1.0
    nx = 256
    ng = 2
    tmax = 0.1
    CFL = 0.8
    # Build the grid
    g = Grid1d_Euler( nx , ng , bc='outflow' )
    # physical indices
    physical = g.physical
    # simulation object
    s = Simulation( g )
    # Assign initial conditions
    initial_condition = "rarefaction"
    s.set_ICs( initial_condition )
    # Save initial conditions
    uinit = s.grid.U.copy()
    # evolve forward in time
    s.evolve( CFL , tmax )
    # get the final values of variables that live on grid
    g = s.grid
    # plot them
    fig,ax = plt.subplots(3,1,sharex=True)
    for var in range(g.NVAR):
        ax[var].plot( g.xs[physical],\
                g.U[var,physical], color='k',\
                label='t={}'.format(tmax))
        ax[var].plot( g.xs[physical] ,uinit[var,physical],\
                ls=":",color="red",zorder=-1,\
                label='initial configuration')
        ax[var].legend(loc='best')
        ax[var].set_ylabel(g.cons_vars[var])

    plt.xlabel("$x$")
    # plt.ylabel("$u$")
    plt.savefig("plots/fv-burger-{}.pdf".format( initial_condition ) )
    plt.suptitle("Solution for {} wave".format(initial_condition))
    plt.show()


if __name__=="__main__":
    import matplotlib.pyplot as plt
    from grid_1d import Grid1d_Euler
    from simulation import Simulation
    main()