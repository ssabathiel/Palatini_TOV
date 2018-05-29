#include "header/glob_variables.h"

void configure()
{
    GR_theory =         0;
    fR_theory =         0;
    fRQ_theory =        0;
    fR_lim_theory =     0;
    fR_metric_theory =  1;

    analytical_EOS =    0;
    tabular_EOS =       0;
    ply_EOS =           1;

    double alpha = 2*pow(10,9);
    double beta = 1.0*pow(10,-7);

    Rp = 1.0/alpha;
    Rq = 1.0/beta;

}
