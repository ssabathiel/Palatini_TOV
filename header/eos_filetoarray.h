#ifndef EOS_FILETOARRAY
#define EOS_FILETOARRAY
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



void readCSV(std::istream &input, std::vector< std::vector<std::string> > &output);
std::vector< std::vector<std::string> > csv_to_array(const char* i_file_name, const char* o_file_name);
alglib::spline1dinterpolant string_to_double_EOS_por(std::vector< std::vector<std::string> > csvData_main,std::vector<double> &rhos, std::vector<double> &press);
alglib::spline1dinterpolant string_to_double_EOS_rop(std::vector< std::vector<std::string> > csvData_main,std::vector<double> &rhos, std::vector<double> &press);
#endif // EOS_FILETOARRAY

