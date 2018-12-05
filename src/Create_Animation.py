
def initialize_animation():
    line.set_data([], [])
    return line,

def animate(t):
    e = EulerSolver(Nx=400 , a=0.0 , b=1.0 , cfl=0.3, time_order=1,spatial_order=1 )
    e.setSod()
    e.evolve(t)
    line.set_data( e.x,e.W[2,:] )
    return line,

if __name__=="__main__":
    fig = plt.figure()
    ax = plt.axes(xlim=(0, 1), ylim=(-2, 2))
    line, = ax.plot([], [], lw=2)
    anim = animation.FuncAnimation(fig, animate, init_func=initialize_animation,
                                   frames=200, interval=20, blit=True)
    anim.save('Low-Sod-Shock-Tube.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
