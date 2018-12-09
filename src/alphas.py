from numpy import maximum
from eigenvalues import lambdaP,lambdaM
def alphaP( WL , csL , WR , csR ):
    ap = maximum( 0 , lambdaP( WL[1,:] , csL ) )
    ap = maximum( ap , lambdaP( WR[1,:] , csR ) )
    return ap

def alphaM( WL , csL , WR , csR ):
    am = maximum( 0 , -lambdaM( WL[1,:] , csL ) )
    am = maximum( am , -lambdaM( WR[1,:] , csR ) )
    return am
