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

    alpha =1.0*pow(10,-9);
    beta = 10.0*pow(10,9);

    Rp = 1.0/alpha;
    Rq = 1.0/beta;


    alpha_iterate = 0;
    alpha_start =alpha;
    alpha_end   =+12.0*pow(10,+9);
    alpha_step = 0.5*pow(10,+9);

}
