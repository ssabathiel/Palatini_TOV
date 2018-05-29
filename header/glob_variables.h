#ifndef GLOB_VARIABLES
#define GLOB_VARIABLES
#include <vector>
#include <cmath>
#include "linalg_library/interpolation.h"
#include "linalg_library/spline.h"
#include "meta_functions.h"



    extern double Rp;
    extern double Rq;

    extern bool GR_theory;
    extern bool fR_theory;
    extern bool fRQ_theory;
    extern bool fR_lim_theory;
    extern bool fR_metric_theory;

    extern bool analytical_EOS;
    extern bool tabular_EOS;
    extern bool ply_EOS;

    extern int th;
    extern int num_an;

    extern double p_m;

    extern double min_press;



    namespace constants{
        double const pi = 3.141592653589793238462643;   //half a Tau
        long double const ggrav = 1.0;// 6.67408*pow(10,-8);   //ORIGINALLY: -8       //Grav-Konstante: 6.67408 × 10-11 m3 kg-1 s-2
        double const G=ggrav;
        double const clight = 1.0;//2.99792458*pow(10,10);    //Speed of light: 2,99*10^8 m/

        long double const g_grav = 6.67408*pow(10,-8);   //ORIGINALLY: -8       //Grav-Konstante: 6.67408 × 10-11 m3 kg-1 s-2
        double const c_light = 2.99792458*pow(10,10);
        double const Gdc2 = g_grav/(c_light*c_light);
        double const Gdc4 = g_grav/(c_light*c_light*c_light*c_light);

        double const g_c2 = ggrav/pow(clight,2);
        double const msun = 1.98892*pow(10,33);         //Mass of sun: 1.989*10^30 kg
        double const rsun = 7*pow(10,10);
        //EOS-constants
        double const h_bar = 1.05457*pow(10,-34)*pow(10,7);                                                          //or: 10,7
        double const m_n = 1.674927*pow(10,-24);              //Mass of neutron: 1.67*10^-27 kg                       or: 10,-24

        //Range and accuracy
        double const delta = 0.23;                      // Scalar for adapted integration-step dr   //original 0.23
        double const rho_center = 1.0*pow(10,14);       // rho_center and value if rho=const
        double const rho_center_max2 = 6.749*pow(10,15); //pow(10,20);
        double const max_r_iterations=10000;
        //double const min_press = pow(1*10,22);  //10,-8         // Set minimal Pressure, otherwise numerical issues, -9 or +26
        double const r_0 = pow(10,-10);
        double const a =  1.0;
        double const kappa_2 = 8*pi;//8*pi*ggrav/pow(clight,2);
        double const gamma_0 = 5.0/3.0;
        double const gamma_1 = 6.0/2.0;            //or: 2.5 or 3.0
        double const rho_1 = pow(10,17);            // or: 10,17
        double const K_0 = 0.2*pow(3*pi*pi,2.0/3.0)*h_bar*h_bar/pow(m_n,8.0/3.0);     //pow(10,8)
        double const gammma = 1.7;                       // For exponential EOS
        double const K = 5.;
        double const p_1 = K_0*pow(rho_1,gamma_0);
        double const K_1 = p_1/pow(rho_1, gamma_1);

    }




    extern int ccount;
    extern int r_count;

    extern alglib::spline1dinterpolant rho_of_press_alg;
    extern alglib::spline1dinterpolant press_of_rho_alg;

    extern tk::spline pf;       //pf = p(rho)-function
    extern tk::spline rhof;     //rhof = rho(p)-function

    extern std::vector<FILE*> p_rho_profiles;


    extern std::vector<rho_of_p_functions> rho_of_p;
    extern std::vector<p_of_rho_functions> p_of_rho;
    extern std::vector<drho_dp_functions> drho_dp;
    extern std::vector<ddrho_dpp_functions> ddrho_dPP;








#endif // GLOB_VARIABLES

