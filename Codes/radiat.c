#include "pluto.h"
#include <math.h>
#include <stdio.h>

#define THOMSON_CROSS_SECTION 6.65e-25
#define LIGHT_SPEED 2.99792458e10
#define FINE_STRUCTURE_CONSTANT 7.2973525698e-3
#define BOLTZMANN_CONSTANT 1.3806488e-16
#define ELECTRON_MASS 9.10938291e-28
#define HIGH_T_ELECTRON_FRACTION 1.1792


double CoolingFunction(double T, double nH, double Z)
{
    double n=log10(nH);

    double lh_a = 4.86567e-13;
    double lh_b = -2.21974;
    double lh_c = 1.35332e-5;
    double lh_d = 9.64775;
    double lh_e = 1.11401e-9;
    double lh_f = -2.66528;
    double lh_g = 6.91908e-21;
    double lh_h = -0.571255;
    double lh_i = 2.45595e-27;
    double lh_j = 0.49521;

    double Lambda_h = (lh_a * pow(T,lh_b) + pow(lh_c * T,lh_d)*(lh_e * pow(T,lh_f) + lh_g * pow(T,lh_h))) / (1 + pow(lh_c * T, lh_d)) + lh_i * pow(T,lh_j);

    double dh_a = 2.84738;
    double dh_b = 3.62655e13;
    double dh_g1 = -3.18564e-4;
    double dh_g2 = 4.8323e-3;
    double dh_g3 = -0.0225974;
    double dh_g4 = 0.0245446;

    double dh_g = dh_g1 * pow(n,4) + dh_g2 * pow(n,3) + dh_g3 * pow(n,2) + dh_g4 * n ;
	
    double D_h = ( pow(T,dh_a) + dh_b * dh_g )/( pow(T,dh_a) + dh_b );
    
    // eq 8
    double mh_a = -0.000847633;
    double mh_b = 0.0127998;
    double mh_c = 45209.3;
    double mh_d = 2.92145e8;

    double M_h=(Z-1.)*(mh_a*n+mh_b)*exp(-(pow((T-mh_c),2.))/mh_d)+1.; 

    // eq 9
    double lm_a = 6.88502e30;
    double lm_b = -1.90262;
    double lm_c = 2.48881e17;
    double lm_d = 0.771176;
    double lm_e = 3.00028e-28;
    double lm_f = 0.472682;

    double Lambda_m = (lm_e * pow(T,lm_f) + pow((lm_a * pow(T, lm_b) + lm_c * pow(T, lm_d)), -1));

    //Eq 10 of Wang et al.2014
    double dm_a = 3.29383;
    double dm_b = 8.82636e14;
    double dm_g1 = 0.00221438;
    double dm_g2 = -0.0353337;
    double dm_g3 = 0.0524811;
	
    double D_m = ( pow(T,dm_a) + dm_b * (dm_g1 * pow(n,3) + dm_g2 * pow(n,2) + dm_g3 * n + 1) )/( pow(T,dm_a) + dm_b );

    // Eq 13 of Wang et al.2014
    double Lambda_e = ((HIGH_T_ELECTRON_FRACTION * FINE_STRUCTURE_CONSTANT * THOMSON_CROSS_SECTION * pow(BOLTZMANN_CONSTANT,2))/(ELECTRON_MASS * LIGHT_SPEED)) * 2.63323e3 * pow(T,1.708064);
	
    double me_a = 0.00769985;
    double me_b = 24683.1;
    double me_c = 0.805234;

    // eq 14
    double M_e = ((me_a*Z-me_a+1)*(pow(T,me_c))+me_b)/(pow(T,me_c)+me_b);

    double Lambda = D_h * Lambda_h * M_h + D_m * Z * Lambda_m + M_e * Lambda_e;
    return Lambda;
}

void Radiat (double *v, double *rhs)
{
  g_minCoolingTemp=2.e4;
  double  mu, T, Tmid, scrh, dT, prs;
  static double E_cost;

  E_cost = UNIT_LENGTH/UNIT_DENSITY/pow(UNIT_VELOCITY, 3.0);

  double nH=v[RHO]/1.e3;
  double Z=1.;
  prs = v[RHOE]*(g_gamma-1.0);
  if (prs < 0.0) {
    prs     = g_smallPressure;
    v[RHOE] = prs/(g_gamma - 1.0);
  }

  mu  = MeanMolecularWeight(v);
  T   = prs/v[RHO]*KELVIN*mu;
  printf("%g %g %g\n",T,nH,prs);

  scrh     = CoolingFunction(T,nH,Z);
  rhs[RHOE] = -scrh*v[RHO]*v[RHO];
    
  scrh       = UNIT_DENSITY/(CONST_amu*mu);  
  rhs[RHOE] *= E_cost*scrh*scrh;
}


/* ***************************************************************** */
// void Radiat (double *v, double *rhs)
/*!
 *   Provide r.h.s. for tabulated cooling.
 * 
 ******************************************************************* */
// {
//   int    klo, khi, kmid;
//   static int ntab;
//   double  mu, T, Tmid, scrh, dT, prs;
//   static double *L_tab, *T_tab, E_cost;
  
//   FILE *fcool;

// /* -------------------------------------------
//         Read tabulated cooling function
//    ------------------------------------------- */

//   if (T_tab == NULL){
//     print (" > Reading table from disk...\n");
//     fcool = fopen("cooltable.dat","r");
//     if (fcool == NULL){
//       print ("! Radiat: cooltable.dat could not be found.\n");
//       QUIT_PLUTO(1);
//     }
//     L_tab = ARRAY_1D(20000, double);
//     T_tab = ARRAY_1D(20000, double);

//     ntab = 0;
//     while (fscanf(fcool, "%lf  %lf\n", T_tab + ntab, 
//                                        L_tab + ntab)!=EOF) {
//       ntab++;
//     }
//     E_cost = UNIT_LENGTH/UNIT_DENSITY/pow(UNIT_VELOCITY, 3.0);
//   }

// /* ---------------------------------------------
//             Get pressure and temperature 
//    --------------------------------------------- */

//   prs = v[RHOE]*(g_gamma-1.0);
//   if (prs < 0.0) {
//     prs     = g_smallPressure;
//     v[RHOE] = prs/(g_gamma - 1.0);
//   }

//   mu  = MeanMolecularWeight(v);
//   T   = prs/v[RHO]*KELVIN*mu;

//   if (T != T){
//     printf (" ! Nan found in radiat \n");
//     printf (" ! rho = %12.6e, prs = %12.6e\n",v[RHO], prs);
//     QUIT_PLUTO(1);
//   }

//   if (T < g_minCoolingTemp) { 
//     rhs[RHOE] = 0.0;
//     return;
//   }

// /* ----------------------------------------------
//         Table lookup by binary search  
//    ---------------------------------------------- */

//   klo = 0;
//   khi = ntab - 1;

//   if (T > T_tab[khi] || T < T_tab[klo]){
//     print (" ! T out of range   %12.6e\n",T);
//     QUIT_PLUTO(1);
//   }

//   while (klo != (khi - 1)){
//     kmid = (klo + khi)/2;
//     Tmid = T_tab[kmid];
//     if (T <= Tmid){
//       khi = kmid;
//     }else if (T > Tmid){
//       klo = kmid;
//     }
//   }

// /* -----------------------------------------------
//     Compute r.h.s
//    ----------------------------------------------- */

//   dT       = T_tab[khi] - T_tab[klo];
//   scrh     = L_tab[klo]*(T_tab[khi] - T)/dT + L_tab[khi]*(T - T_tab[klo])/dT;
//   rhs[RHOE] = -scrh*v[RHO]*v[RHO];
  
//   scrh       = UNIT_DENSITY/(CONST_amu*mu);  
//   rhs[RHOE] *= E_cost*scrh*scrh;
// }