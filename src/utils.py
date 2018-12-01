def minmod( x , y , z ):
    return( 1./4. * np.fabs( np.sign(x) + np.sign(y)) * \
            (np.sign(x) + np.sign(z)) * \
            min(np.minimum(np.minimum(np.fabs(x),np.fabs(y)),np.fabs(z))))

def compute_l1_error( numerical , exact , deltaX ):
    diff = np.abs( numerical - exact )
    summ = diff.sum()
    return deltaX * summ

def fit_line( xdata , ydata ):
    # fit data to a line
    b, m = polyfit( xdata , ydata , 1 )
    return b , m

def line( x , b , m ):
    return m * x + b

def plot_convergence(order='low'):
    """
    this function plots the convergence rate of then
    scheme
    """
    Ns = [ 32 , 64 , 128 , 256 , 512 ]
    t=0.1
    x0=0.5
    a=0.0
    b=1.0
    rhol = 1.0; vl = 0.0; pl = 1.0
    rhor = 0.1; vr = 0.0; pr = 0.125
    gamma=1.4

    rerr = []
    verr = []
    perr = []

    L1 = np.zeros( len(Ns) )
    for i in range(len(Ns)):
        deltaX = (b-a)/Ns[i]
        if order=='low':
            e = EulerSolver(Nx = Ns[i],time_order=1,spatial_order=1)
        else:
            e = EulerSolver(Nx = Ns[i],time_order=2,spatial_order=2)
        e.setSod()
        e.evolve(t)
        xexact , rexact , vexact , pexact = riemann( a , b , x0 , Ns[i] , t ,rhol ,vl , pl , rhor , vr , pr , gamma )
        rerr.append( compute_l1_error( e.W[0,:] , rexact , deltaX ))
        verr.append( compute_l1_error( e.W[1,:] , vexact , deltaX ))
        perr.append( compute_l1_error( e.W[2,:] ,   pexact , deltaX ))

    print("L1 Errors on Density",rerr)
    print("L1 Errors on Velocity",verr)
    print("L1 Errors on Pressure",perr)

    logns = np.log( Ns )
    logrs = np.log( rerr )
    logvs = np.log( verr )
    logps = np.log( perr )

    rb , rm = fit_line( logns , logrs )
    vb , vm = fit_line( logns , logvs )
    pb , pm = fit_line( logns , logps )

    print("Density slope:" , rm)
    print("Velocity slope:" , vm)
    print("Pressure slope:" , pm)

    fig = plt.figure()

    plt.plot( logns, logrs , color='red' , marker='o',label="L1 Error on Density")
    plt.plot( logns, logvs , color='green' , marker='o',label="L1 Error on Velocity")
    plt.plot( logns, logps , color='blue' , marker='o',label="L1 Error on Pressure")

    plt.plot( logns , line(logns , rb , rm ) , linestyle='dashed', color='red' , label='L1 Error Fit on Density')
    plt.plot( logns , line(logns , vb , vm ) , linestyle='dashed',color='green' , label='L1 Error Fit on Velocity')
    plt.plot( logns , line(logns , pb , pm ) , linestyle='dashed',color='blue' , label='L1 Error Fit on Pressure')

    plt.title("Convergence Plot for Sod Shock Tube")
    plt.xlabel("$log(N)$")
    plt.ylabel("$log(L_1)$")
    plt.legend()

def f(x,x0,sigma):
    return(np.where(abs(x-x0)<sigma , (1-((x-x0)/sigma)**2)**2 , 0))
