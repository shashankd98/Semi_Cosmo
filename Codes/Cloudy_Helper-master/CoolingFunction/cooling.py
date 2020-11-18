#! /usr/bin/env python3

import sys, math

def CoolingFunction(T, n, Z):

    # Thomson cross-section, cm^2
    THOMSON_CROSS_SECTION = 6.65e-25

    # light speed, cm/s
    LIGHT_SPEED = 2.99792458e10

    # fine structure constant
    FINE_STRUCTURE_CONSTANT = 7.2973525698e-3

    # boltzmann constant, erg/K
    BOLTZMANN_CONSTANT = 1.3806488e-16

    # electron mass, gram
    ELECTRON_MASS = 9.10938291e-28

    # electron density/hydrogen density at solar abudance and
    # fully ionized status, value here derived with abundance
    # GASS discribed in CLOUDY */
    HIGH_T_ELECTRON_FRACTION = 1.1792


    # Eq 6 of Wang et al.2014 
    lh_a = 4.86567e-13
    lh_b = -2.21974
    lh_c = 1.35332e-5
    lh_d = 9.64775
    lh_e = 1.11401e-9
    lh_f = -2.66528
    lh_g = 6.91908e-21
    lh_h = -0.571255
    lh_i = 2.45595e-27
    lh_j = 0.49521

    Lambda_h = (lh_a * math.pow(T,lh_b) + math.pow(lh_c * T,lh_d)*(lh_e * math.pow(T,lh_f) + lh_g * math.pow(T,lh_h))) / (1 + math.pow(lh_c * T, lh_d)) + lh_i * math.pow(T,lh_j)

    # Eq 7 of Wang et al.2014 
    dh_a = 2.84738
    dh_b = 3.62655e13
    dh_g1 = -3.18564e-4
    dh_g2 = 4.8323e-3
    dh_g3 = -0.0225974
    dh_g4 = 0.0245446

    dh_g = dh_g1 * math.pow(n,4) + dh_g2 * math.pow(n,3) + dh_g3 * math.pow(n,2) + dh_g4 * n + 1
	
    D_h = ( math.pow(T,dh_a) + dh_b * dh_g )/( math.pow(T,dh_a) + dh_b )

    # Eq 9 of Wang et al.2014 
    lm_a = 6.88502e30
    lm_b = -1.90262
    lm_c = 2.48881e17
    lm_d = 0.771176
    lm_e = 3.00028e-28
    lm_f = 0.472682

    Lambda_m = (lm_e * math.pow(T,lm_f) + math.pow((lm_a * math.pow(T, lm_b) + lm_c * math.pow(T, lm_d)), -1))

    # Eq 10 of Wang et al.2014
    dm_a = 3.29383
    dm_b = 8.82636e14
    dm_g1 = 0.00221438
    dm_g2 = -0.0353337
    dm_g3 = 0.0524811
	
    D_m = ( math.pow(T,dm_a) + dm_b * (dm_g1 * math.pow(n,3) + dm_g2 * math.pow(n,2) + dm_g3 * n + 1) )/( math.pow(T,dm_a) + dm_b )

    # Eq 13 of Wang et al.2014
    Lambda_e = ((HIGH_T_ELECTRON_FRACTION * FINE_STRUCTURE_CONSTANT * THOMSON_CROSS_SECTION * math.pow(BOLTZMANN_CONSTANT,2))/(ELECTRON_MASS * LIGHT_SPEED)) * 2.63323e3 * math.pow(T,1.708064)
	
    me_a = 0.00769985

    # High Temperature Approximation of Eq 14 of Wang et al.2014 
    M_e = me_a * (Z - 1) + 1

    Lambda = D_h * Lambda_h +  D_m * Z * Lambda_m + M_e * Lambda_e

    return Lambda


def CoolingRate(T, n, Z):
    Lambda = CoolingFunction(T, n, Z)
    
    # Election fraction, Eq 12 of Wang et al.2014
    E = 2.1792 - math.exp(3966.27/T)
    a = 0.00769985
    b = 24683.1
    c = 0.805234

    # Eq 14 of Wang et al.2014
    M = ((a * Z - a + 1) * math.pow(T,c) + b)/(math.pow(T,c) + b)

    rate = E * M * math.pow(math.pow(10,n),2) * Lambda
    return rate

def main():
    if( len(sys.argv) != 4 ):
        sys.exit('ERROR: I don\'t get all argument I need!\nPlease provide temperature, density, and metallicity!\n')

    T = float(sys.argv[1])
    if( T < 1e4 or T > 1e10 ):
        print('Caution: Temperature is out of fit range!\n')

    n = float(sys.argv[2])
    if( n > 10 ):
        print('Caution: Density is out of fit range!\n')

    Z = float(sys.argv[3])
    if( Z > 30 ):
        print('Caution: Metallicity is out of fit range!\n')

    CF = CoolingFunction(T, n, Z)
    CR = CoolingRate(T, n, Z)

    print('CoolingFunction is {0:e} and CoolingRate is {1:e}'.format(CF, CR))


if __name__ == '__main__':
    main()
