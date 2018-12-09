from numpy import zeros
def HLLE_Flux( UL, UR , FL , FR , am , ap ):
    FHLL = zeros((3,ap.shape[0]))
    Flux_Difference = zeros((3,ap.shape[0]-1))
    FHLL[:,:] = ( ap * FL + am * FR - ap * am * ( UR - UL ) ) / ( ap + am )
    Flux_Difference = -( FHLL[:,1:] - FHLL[:,:-1] )
    return Flux_Difference
