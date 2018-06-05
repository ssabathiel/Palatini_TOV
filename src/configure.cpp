#include "header/glob_variables.h"

void configure()
{
    GR_theory =         0;
    fR_theory =         0;
    fRQ_theory =        1;
    fR_lim_theory =     0;
    fR_metric_theory =  0;

    analytical_EOS =    1;
    tabular_EOS =       0;
    ply_EOS =           0;
    fps_EOS =           0;
    ap4_EOS =           0;

    double alpha =2.0*pow(10,9);
    double beta = 1.0*pow(10,-20);

    Rp = 1.0/alpha;
    Rq = 1.0/beta;

}
