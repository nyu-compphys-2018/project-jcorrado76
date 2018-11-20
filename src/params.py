
# sod problem
sod_params = {\
        'rho_l' : 1.0 ,\
        'u_l'   : 0.0 ,\
        'p_l'   : 1.0 ,\
        'rho_r' : 0.125,\
        'u_r'   : 0.0  ,\
        'p_r'   : 0.1 ,\
        'tmax'  : 0.2 ,\
        'gamma' : 1.4 ,\
        'cfl'   : 0.8 ,\
        }

# double rarefaction
params = {
    'rho_l': 1.0,
    'u_l': -2.0,
    'p_l': 0.4,
    'rho_r': 1.0,
    'u_r': 2.0,
    'p_r': 0.4,
    'tmax': 0.15,
    'gamma': 1.4,
    'cfl': 0.8
    }

# strong shock
params = {
    'rho_l': 1.0,
    'u_l': 0.0,
    'p_l': 1000.0,
    'rho_r': 1.0,
    'u_r': 0.0,
    'p_r': 0.01,
    'tmax': 0.012,
    'gamma': 1.4,
    'cfl': 0.8
    }

# slow moving shock
params = {
    'rho_l': 5.6698,
    'u_l': -1.4701,
    'p_l': 100.0,
    'rho_r': 1.0,
    'u_r': -10.5,
    'p_r': 1.0,
    'tmax': 1.0,
    'gamma': 1.4,
    'cfl': 0.8
    }
