#include "cooling.h"

/* Calculate Cooling Function */
float Cooling::CoolingFunction(float T, float n, float Z)
{
	using namespace Cooling;
	
	float Lambda;

	/* Eq 6 of Wang et al.2014 */
	float lh_a = 4.86567e-13;
	float lh_b = -2.21974;
	float lh_c = 1.35332e-5;
	float lh_d = 9.64775;
	float lh_e = 1.11401e-9;
	float lh_f = -2.66528;
	float lh_g = 6.91908e-21;
	float lh_h = -0.571255;
	float lh_i = 2.45595e-27;
	float lh_j = 0.49521;

	float Lambda_h = (lh_a * pow(T,lh_b) + pow(lh_c * T,lh_d)*(lh_e * pow(T,lh_f) + lh_g * pow(T,lh_h))) / (1 + pow(lh_c * T, lh_d)) + lh_i * pow(T,lh_j);

	//printf("Lambda_h is %e\n", Lambda_h);

	/* Eq 7 of Wang et al.2014 */
	float dh_a = 2.84738;
	float dh_b = 3.62655e13;
	float dh_g4 = -3.18564e-4;
	float dh_g3 = 4.8323e-3;
	float dh_g2 = -0.0225974;
	float dh_g1 = 0.0245446;

	float dh_g = dh_g4 * pow(n,4) + dh_g3 * pow(n,3) + dh_g2 * pow(n,2) + dh_g1 * n + 1;
	
	float D_h = ( pow(T,dh_a) + dh_b * dh_g )/( pow(T,dh_a) + dh_b );

	//printf("D_h is %e\n", D_h);

	/* Eq 9 of Wang et al.2014 */
	float lm_a = 6.88502e30;
	float lm_b = -1.90262;
	float lm_c = 2.48881e17;
	float lm_d = 0.771176;
	float lm_e = 3.00028e-28;
	float lm_f = 0.472682;

	float Lambda_m = Z * (lm_e * pow(T,lm_f) + pow((lm_a * pow(T, lm_b) + lm_c * pow(T, lm_d)), -1));

	//printf("Lambda_m is %e\n", Lambda_m);

	/* Eq 10 of Wang et al.2014 */
	float dm_a = 3.29383;
	float dm_b = 8.82636e14;
	float dm_g3 = 0.00221438;
	float dm_g2 = -0.0353337;
	float dm_g1 = 0.0524811;
	
	float D_m = ( pow(T,dm_a) + dm_b * (dm_g3 * pow(n,3) + dm_g2 * pow(n,2) + dm_g1 * n + 1) )/( pow(T,dm_a) + dm_b );

	//printf("D_m is %e\n", D_m);

	/* Eq 11 of Wang et al.2014 */
	float Lambda_e = ((HIGH_T_ELECTRON_FRACTION * FINE_STRUCTURE_CONSTANT * THOMSON_CROSS_SECTION * pow(BOLTZMANN_CONSTANT,2))/(ELECTRON_MASS * LIGHT_SPEED)) * 2.63323 * pow(T,1.708064);
	
	float me_a = 0.00769985;

	/* High Temperature Approximation of Eq 14 of Wang et al.2014 */
	float M_e = me_a * (Z - 1) + 1;

	Lambda = D_h * Lambda_h +  Z * D_m * Lambda_m + M_e * Lambda_e;

	return Lambda;
}

float Cooling::CoolingRate(float T, float n, float Z)
{
	using namespace Cooling;
	
	float rate;
	float Lambda = CoolingFunction(T, n, Z);

	/* Election fraction, Eq 12 of Wang et al.2014 */
	float E = 2.1792 - exp(3966.27/T);
	float a = 0.00769985;
	float b = 24683.1;
	float c = 0.805234;

	/* Eq 14 of Wang et al.2014 */
	float M = ((a * Z - a + 1) * pow(T,c) + b)/(pow(T,c) + b);


	rate = E * M * pow(pow(10,n),2) * Lambda;

	return rate;
}
