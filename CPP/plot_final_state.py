if __name__=="__main__":
    import sys 
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np

    try:
        inFileName = sys.argv[1]
    except:
        inFileName = "outputFile.csv"
    data =\
    pd.read_csv( inFileName ,sep=',',header=None,names=['x',r'$\rho$',r'v',r'p'],\
            dtype=np.float64)

    data.plot(x='x',subplots=True, sharex=True, title="Sod Shock Tube to time t=0.1",
            grid=True,legend=True)

    plt.show()
    plt.savefig("sod_shock_tube.jpg",format='jpg',dpi=400,bbox_inches='tight')
