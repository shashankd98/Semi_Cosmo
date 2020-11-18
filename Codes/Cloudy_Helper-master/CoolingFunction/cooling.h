#ifndef COOLING_H_
#define COOLING_H_

#include "stdio.h"
#include "math.h"

/* to avoid name corruption, all variables and funcitons are
 * using namespace 'Cooling' now */
namespace Cooling
{
	/* Replace these constants with yours, if you have difined
	 * them in your code */
	/* Thomson cross-section, cm^2 */
	const float THOMSON_CROSS_SECTION = 6.65e-25;
	
	/* light speed, cm/s */
	const float LIGHT_SPEED = 2.99792458e10;
	/* fine structure constant */
	const float FINE_STRUCTURE_CONSTANT = 7.2973525698e-3;

	/* boltzmann constant, erg/K */
	const float BOLTZMANN_CONSTANT = 1.3806488e-16;

	/* electron mass, gram */
	const float ELECTRON_MASS = 9.10938291e-28;
	/* electron density/hydrogen density at solar abudance and
	 * fully ionized status, value here derived with abundance
	 * GASS discribed in CLOUDY */
	const float HIGH_T_ELECTRON_FRACTION = 1.1792;

	/* change type to what you use, T is temperature (K), n is log(n_H),
	 * n_H in unit of cm-3, Z is metallicity in unit of Z_sun
	 * CoolingFunction return a cooing function in unit of erg cm3 s-1*/
	float CoolingFunction(float T, float n, float Z);

	/* CoolingRate return cooling rate in unit of erg cm-3 s-1 */
	float CoolingRate(float T, float n, float Z);
}
#endif
