#include "header/glob_variables.h"

void configure()
{
    GR_theory =         0;
    fR_theory =         0;
    fRQ_theory =        1;
    fR_lim_theory =     0;
    fR_metric_theory =  0;

    analytical_EOS =    0;
    tabular_EOS =       0;
    ply_EOS =           0;
    fps_EOS =           0;
    ap4_EOS =           1;

    alpha =1.0*pow(10,-9);
    beta = 0.1*pow(10,+9);

    Rp = 1.0/alpha;
    Rq = 1.0/beta;


    alpha_iterate = 0;
    alpha_start =alpha;
    alpha_end   =+2.0*pow(10,+9);
    alpha_step = 0.1*pow(10,+9);


    beta_iterate = 1;
    beta_start =beta;
    beta_end   =+8*pow(10,+9);
    beta_step = 0.2*pow(10,+9);

}
