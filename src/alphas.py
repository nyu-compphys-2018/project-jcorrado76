def alphaP( WL , csL , WR , csR ):
    ap = np.maximum( 0 , self.lambdaP( WL[1,:] , csL ) )
    ap = np.maximum( ap , self.lambdaP( WR[1,:] , csR ) )
    return ap

def alphaM( WL , csL , WR , csR ):
    am = np.maximum( 0 , -self.lambdaM( WL[1,:] , csL ) )
    am = np.maximum( am , -self.lambdaM( WR[1,:] , csR ) )
    return am
