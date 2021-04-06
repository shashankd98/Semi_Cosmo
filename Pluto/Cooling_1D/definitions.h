#define  PHYSICS                        HD
#define  DIMENSIONS                     1
#define  GEOMETRY                       SPHERICAL
#define  BODY_FORCE                     VECTOR
#define  COOLING                        TABULATED
#define  RECONSTRUCTION                 PARABOLIC
#define  TIME_STEPPING                  RK3
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            0

/* -- physics dependent declarations -- */

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */


/* [Beg] user-defined constants (do not change this line) */

#define  SHOCK_FLATTENING               ONED
#define  CHAR_LIMITING                  YES
#define  INTERNAL_BOUNDARY              YES
#define  UNIT_DENSITY                   1.67262171e-27
#define  UNIT_LENGTH                    3.0856775807e21
#define  UNIT_VELOCITY                  1.e7

/* [End] user-defined constants (do not change this line) */
