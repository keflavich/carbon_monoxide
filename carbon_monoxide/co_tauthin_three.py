"""
Self-contained version, copied from
https://github.com/keflavich/w51_singledish_h2co_maps/blob/676ed941272bd41a6ff900ec3f80dbcec44bef69/common_constants.py
"""
# CO constants
coprops12CO10 = {'Aul':7.203e-8*u.Hz, 'Be': 57.63596828e9*u.Hz, 'nu': 115.271202e9 *u.Hz, 'Ju': 1, 'h2toco':1e4}
coprops13CO10 = {'Aul':6.294e-8*u.Hz, 'Be': 55.1010138e9 *u.Hz, 'nu': 110.201370e9 *u.Hz, 'Ju': 1, 'h2toco':6e5}
coprops12CO21 = {'Aul':6.910e-7*u.Hz, 'Be': 57.63596828e9*u.Hz, 'nu': 230.5380e9   *u.Hz, 'Ju': 2, 'h2toco':1e4}
coprops13CO21 = {'Aul':6.038e-7*u.Hz, 'Be': 55.1010138e9 *u.Hz, 'nu': 220.3986765e9*u.Hz, 'Ju': 2, 'h2toco':6e5}
coprops12CO32 = {'Aul':2.497e-6*u.Hz, 'Be': 57.63596828e9*u.Hz, 'nu': 345.7959899e9*u.Hz, 'Ju': 3, 'h2toco':1e4}
coprops13CO32 = {'Aul':2.181e-6*u.Hz, 'Be': 55.1010138e9 *u.Hz, 'nu': 330.5879601e9*u.Hz, 'Ju': 3, 'h2toco':6e5}


def nj(tem,Ju,numlevs=50,Be=coprops13CO10['Be']):
    """
    Energy level population of upper state for CO
    """
    from astropy.constants import h, k_B
    J = np.arange(numlevs)
    ntotDn0 = np.sum( (2*J+1)*np.exp(-J*(J+1.0)*Be*h/(k_B*tem)) )
    if Ju == 0:
        return 1.0/ntotDn0
    expBe = np.exp(-h*Be*(Ju*(Ju+1))/ ( k_B * tem ) )
    nuDn0 = (2*Ju+1.0) * expBe
    #ntotDnu = ntotDn0 * (2.0*Jl+1.0)/(2.0*Ju+1.0) * exp(h*Be*(Ju*(Ju+1))/ ( k_B * tem ) ) 
    #ntotDnu = ntotDn0 * ((2.0*Jl+1.0)/(2.0*Ju+1.0) * exp(h*Be*(Ju*(Ju+1))/ ( k_B * tem ) ) - \
    #                    (2.0*Jl+1.0)/(2.0*Ju+1.0) * exp(h*Be*(Ju*(Ju+1))/ ( k_B * TCMB*u.K ) ) )
    return nuDn0/ntotDn0

# CO to H2 conversion factor
# Taken from co_tauthin.py (other repository)
def cotocol(tex=20,coprops=coprops13CO10):
    """
    Conversion ratio - ratio of particular CO molecule to H2
    """
    from astropy.constants import h, k_B, c
    from numpy import pi
    nu = coprops['nu']
    Ju = coprops['Ju']
    Be = coprops['Be']
    Aul = coprops['Aul']
    conversionratio = coprops['h2toco']
    if not isinstance(tex,np.ndarray):
        if hasattr(tex,'__len__') and len(tex) > 1:
            tex = np.array(tex)
        else:
            tex = np.array([tex])
    if not hasattr(tex,'unit') or tex.unit is None:
        tex = tex * u.K
    ntotDnu = 1.0 / np.array([ nj(T,Ju,Be=Be,numlevs=50) for T in tex ]) 
    expcmb = np.exp(h*nu/(k_B*TCMB*u.K))
    Nupper = ( 8 * pi * nu**2 ) / ( c**3 * Aul ) * k_B / h * (expcmb - 1) / ( expcmb - (np.exp(h*nu/(k_B*tex))) ) # - (exp(hplanck*nu/(k_B*TCMB*u.K))-1)**-1)
    return (conversionratio * ntotDnu * Nupper).to(u.cm**-2 / (u.K * u.km/u.s))

