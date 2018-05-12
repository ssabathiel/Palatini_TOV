#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>
#include <cmath>
#include <cstdlib>  // for atof string to double conversion
#include "linalg_library/spline.h"
#include <iomanip>
#include "linalg_library/interpolation.h"
#include "header/plotradmassrho.h"
#include "header/eos_filetoarray.h"
#include "header/eos_functions.h"
#include "header/glob_variables.h"
#include "header/get_gradient_gr.h"

using namespace constants;
using namespace std;




pair<double, double> get_gradients_GR(double m, double press, double r)
{
    if(press<min_press) { press = min_press; }
    else                { press = press;     }

    double rho = rho_of_p[num_an](press);

    rho *= Gdc2;
    press *= Gdc4;
    m = m *= Gdc2;

    double dm_dr = 4*pi*rho*pow(r,2);
    double dp_dr;

    //Get dp/dr = ...TOV...
    if(r>1*pow(10,-6)) { dp_dr = (m +4.0*pi*r*r*r*press/pow(clight,2))/(r*(r-2.0*ggrav*m/pow(clight,2))); }
    else { dp_dr = (4.0*pi*r*press/pow(clight,2)); }

    dp_dr = -ggrav*( rho+press/pow(clight,2) )*dp_dr;

    dm_dr /= Gdc2;
    m /= Gdc2;
    dp_dr /= Gdc4;


    return make_pair(dm_dr, dp_dr);
}

