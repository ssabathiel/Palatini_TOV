#ifndef EOS_FUNCTIONS
#define EOS_FUNCTIONS
#include <vector>
#include <math.h>
#include <cmath>
#include <cstdlib>  // for atof string to double conversion
#include "linalg_library/spline.h"
#include <iomanip>
#include "linalg_library/interpolation.h"
#include "glob_variables.h"


double rho_of_p_num(double p);
double p_of_rho_num(double rho);
double drho_dp_num(double p);
void create_spline(char * file_name);
double ddrho_dPP_num(double pressy);

#endif // EOS_FUNCTIONS

