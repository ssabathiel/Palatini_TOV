#ifndef OUTPUT
#define OUTPUT

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>  // for atof string to double conversio
#include "eos_filetoarray.h"
#include "linalg_library/interpolation.h"
#include "linalg_library/spline.h"
#include "glob_variables.h"



void output_to_arrays(const char* i_file_name,  std::vector<double> &rhos, std::vector<double> &masses, std::vector<double> &radii );
std::string DoublePowToString(double Rpp);
std::string create_theory_ID();
std::string create_theory_ID_wo_alpha_beta();

#endif // OUTPUT

