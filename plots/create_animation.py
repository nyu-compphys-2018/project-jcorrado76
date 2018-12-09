

def initialize_animation():
    line.set_data([], [])
    return line,

def animate(t):
    e = EulerSolver(Nx=400 , a=0.0 , b=1.0 , cfl=0.3, time_order=2,spatial_order=2 )
    e.setSod()
    e.evolve(t)
    line.set_data( e.x,e.W[2,:] )
    return line,
