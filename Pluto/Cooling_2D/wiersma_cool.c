#include "pluto.h"
/* ***************************************************************** */

double solver(double);

double WiersmaCool (double T)
/*!
 *   returns Lambda[T] in cgs units.
 * 
 ******************************************************************* */
{
  int    klo, khi, kmid;
  static int ntab;
  double  Tmid, scrh, dT;
  static double *L_tab, *T_tab, *Z_tab, met;
  FILE *fcool;

  double z=solver(g_time*UNIT_LENGTH/UNIT_VELOCITY);
 // if(T>1.e5) met=pow(10.,-1.+z/6.*(-1.));
 // else met=pow(10.,-2.+z/6.*(-3.));
  met=0.01;
/* -------------------------------------------
        Read tabulated cooling function
   ------------------------------------------- */

  if (T_tab == NULL){
    print (" > Reading table from disk...\n");
    fcool = fopen("wiersma_cooltable.dat","r");
    if (fcool == NULL){
      print ("! LamCool: cooltable.dat could not be found.\n");
      QUIT_PLUTO(1);
    }
    L_tab = ARRAY_1D(20000, double);
    T_tab = ARRAY_1D(20000, double);
    Z_tab=ARRAY_1D(20000,double);

    

    ntab = 0;
    while (fscanf(fcool, "%lf  %lf %lf \n", T_tab + ntab,
                                       L_tab + ntab, Z_tab+ntab)!=EOF) {
      ntab++;
    }
  }

/* ---------------------------------------------
            Make sure that T is well-defined 
   --------------------------------------------- */

  if (T != T){
    printf (" ! Nan found in lam_cool \n");
    printf (" ! T = %12.6e\n", T);
    QUIT_PLUTO(1);
  }
  g_minCoolingTemp=2.e4;
  if (T < g_minCoolingTemp) {
    return 0.0;
  }
  // printf("%g \n",g_minCoolingTemp);
/* ----------------------------------------------
        Table lookup by binary search  
   ---------------------------------------------- */

  klo = 0;
  khi = ntab - 1;

  if (T < T_tab[klo]){
    //print (" T<Tlo    %12.6e\n",T);
    //QUIT_PLUTO(1);
    return L_tab[klo]+met*Z_tab[klo];
  }

  if (T > T_tab[khi]){
    //print (" T>Thi   %12.6e\n",T);
    //QUIT_PLUTO(1);
    return L_tab[khi]+met*Z_tab[khi];
  }

  while (klo != (khi - 1)){
    kmid = (klo + khi)/2;
    Tmid = T_tab[kmid];
    if (T <= Tmid){
      khi = kmid;
    }else if (T > Tmid){
      klo = kmid;
    }
  }

/* -----------------------------------------------
    Compute and return Lambda 
   ----------------------------------------------- */

  dT       = T_tab[khi] - T_tab[klo];
  scrh     = (L_tab[klo]*(T_tab[khi] - T)/dT + L_tab[khi]*(T - T_tab[klo])/dT)+met*(Z_tab[klo]*(T_tab[khi] - T)/dT + Z_tab[khi]*(T - T_tab[klo])/dT);
  // printf("%g %g\n",T,scrh);
  return scrh;
  // return 1.e-23;
}
