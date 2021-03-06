//////////////////////////////////////////////
// Code to Solve the TOV-equation, with variable EOS-mode
//      done by Silvester Sabathiel, supervisor: Dr. Hèlios Sanchis-Alepuz
//      project within Masterthesis, start: April 2017
///////////////////////////////////////////////

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iostream>
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
#include "header/output.h"
#include "header/eos_analytical.h"
#include "header/ply_analytical.h"
#include "header/integrate_star.h"
#include "header/meta_functions.h"
#include "header/configure.h"
#include "header/whole.h"

using namespace std;
using namespace constants;



//////////////////////////
// Define global constants:   everything in cm,g,seconds, cgs
//////////////////////////

int ccount=0;
int r_count=0;


//////////////////////////
//Define global data (arrays)
//////////////////////////


alglib::spline1dinterpolant rho_of_press_alg;
alglib::spline1dinterpolant press_of_rho_alg;


string const EOS_mode = "Spline";
//string const EOS_mode = "Const";
//string const EOS_mode = "Polytrope";

bool plot = 1;


vector<FILE*> p_rho_profiles;
pair<double, double> dm_dp_dr;


///////////////////////
//PRINT OR NOT TO PRINT
///////////////////////
bool print_count = 1;
bool print_EOS_check = 0;
bool print_each_integ_step = 0;
bool print_dp_dm = 0;      //Easier to checkk if realistic or not then derivative
bool print_dr = 0;


bool GR_theory;
bool fR_theory;
bool fR_lim_theory;
bool fR_metric_theory;
bool fRQ_theory;

bool analytical_EOS;
bool ply_EOS;
bool tabular_EOS;
bool fps_EOS;
bool ap4_EOS;

double Rp;
double Rq;
double alpha;
double beta;
double alpha_start;
double alpha_end;
double alpha_step;
bool alpha_iterate;

double beta_start;
double beta_end;
double beta_step;
bool beta_iterate;

int th;
int num_an;
int eos_id;
double min_press;

double a1;
 double a2;
 double a3;
double a4;
double a5;
double a6;
double a7;
double a8;
double a9;
double a10;
double a11;
double a12;
double a13;
double a14;
double a15;
double a16;
double a17;
double a18;

vector<rho_of_p_functions> rho_of_p;
vector<p_of_rho_functions> p_of_rho;
vector<drho_dp_functions> drho_dp;
vector<ddrho_dpp_functions> ddrho_dPP;
vector<int> detailed_stars;

int main()
{
    configure();
    rho_of_p.push_back(rho_of_p_analytical);
    rho_of_p.push_back(rho_of_p_num);
    rho_of_p.push_back(rho_of_p_analytical_ply);

    p_of_rho.push_back(p_of_rho_analytical);
    p_of_rho.push_back(p_of_rho_num);
    p_of_rho.push_back(p_of_rho_analytical_ply);

    drho_dp.push_back(drho_dp_analytical);
    drho_dp.push_back(drho_dp_num);
    drho_dp.push_back(drho_dp_analytical_ply);

    ddrho_dPP.push_back(ddrho_dPP_analytical);
    ddrho_dPP.push_back(ddrho_dPP_num);
    ddrho_dPP.push_back(ddrho_dPP_analytical_ply);


    if(analytical_EOS==1)
    {

    // SLY
    a1 =6.22;
     a2 =6.121;
     a3 =0.005925;
     a4 =0.16326;
     a5 =6.48;
     a6 =11.4971;
     a7 =19.105;
     a8 =0.8938;
     a9 =6.54;
     a10 = 11.4950;
     a11 = -22.775;
     a12 = 1.5707;
     a13 = 4.3;
     a14 = 14.08;
     a15 = 27.80;
     a16 = -1.653;
     a17 = 1.50;
     a18 = 14.67;
    }

    if(fps_EOS==1)
    {
    // FPS
     a1 =6.22;
     a2 =6.121;
     a3 =0.006004;
     a4 =0.16345;
     a5 =6.50;
     a6 =11.8440;
     a7 =17.24;
     a8 =1.065;
     a9 =6.54;
     a10 = 11.8421;
     a11 = -22.003;
     a12 = 1.5552;
     a13 = 9.3;
     a14 = 14.19;
     a15 = 23.73;
     a16 = -1.508;
     a17 = 1.79;
     a18 = 15.13;
   }

    if(ap4_EOS==1)
    {
    // AP4
     a1 =4.3290;
     a2 =4.3622;
     a3 =9.1131;
     a4 =-0.475;
     a5 =3.4614;
     a6 =14.8800;
     a7 =21.3141;
     a8 =0.1023;
     a9 =0.0495;
     a10 = 4.9401;
     a11 = 10.2957;

   }



    if(fps_EOS==1 || analytical_EOS==1 || ap4_EOS==1)
    {
        create_points_from_analytical_EOS("./EOS_data/sly_points", p_of_rho_analytical);
        char* test_dat = "EOS_data/sly_points";    //test_data.dat
        create_spline(test_dat);
    }

    if(GR_theory==1){th=0;}
    else if(fR_theory==1){th=1;}
    else if(fRQ_theory==1){th=2;}
    else if(fR_lim_theory==1){th=3;}
    else if(fR_metric_theory==1){th=4;}

    if(analytical_EOS==1){num_an=0;}
    else if(tabular_EOS==1){num_an=1;}
    else if(ply_EOS==1){num_an=2;}
    else if(fps_EOS==1){num_an=0;}
    else if(ap4_EOS==1){num_an=0;}

    string pal_metr_str = "";
    //if(fR_metric_theory==1){pal_metr_str = "_Met";}
    //if(fR_theory==1){pal_metr_str = "_Pal";}

    string alpha_maxmass_file = "./Plotting/Results/";
    alpha_maxmass_file.append(create_theory_ID_wo_alpha_beta() );
    alpha_maxmass_file.append("_Maxmass_alpha");
    FILE *ifp_maxmass_alpha=fopen(alpha_maxmass_file.c_str(),"w");
    rewind(ifp_maxmass_alpha);

    alpha_maxmass_file = "./Plotting/Results/";
    alpha_maxmass_file.append(create_theory_ID_wo_alpha_beta() );
    alpha_maxmass_file.append("_Maxmass_beta");
    FILE *ifp_maxmass_beta=fopen(alpha_maxmass_file.c_str(),"w");
    rewind(ifp_maxmass_beta);

    alpha_maxmass_file = "./Plotting/Results/";
    alpha_maxmass_file.append(create_theory_ID_wo_alpha_beta() );
    alpha_maxmass_file.append("_Maxmass_alpha_beta");
    FILE *ifp_maxmass_alpha_beta=fopen(alpha_maxmass_file.c_str(),"w");
    rewind(ifp_maxmass_alpha_beta);

    FILE *ifp;
    FILE *ifp2;

    //th=0;

    alpha = alpha_start;
    beta = beta_start;
    while(alpha<alpha_end)
    {
        cout << "ALPHA= " << alpha << endl;
        while(beta<beta_end)
        {
            cout << "BETA= " << beta << endl;

            double Mass;
            double Radius;
            vector<double> Masses;
            vector<double> Radi;

            double rho_center_now = rho_center;
            double rho_center_max = rho_center_max2;
            double Maxmass = 0;
            double Maxmass_Radius = 0;

             if(abs(alpha)>pow(10,-80) && abs(beta)>pow(10,-80))
             {
                ///////////////////////////////////////
                //LOOP OVER RHO_CENTERS if multiple stars, other wise INTEGRATE ONLY ONCE
                ///////////////////////////////////////


                string theory_ID = create_theory_ID();
                string output_path = "./Plotting/Results/TOV_output_";
                output_path.append(theory_ID);
                output_path.append(pal_metr_str);
                ifp=fopen(output_path.c_str(),"w");
                rewind(ifp);
                ifp2=fopen("TOV_output","w");
                rewind(ifp2);


                //vector<int> detailed_stars;
                detailed_stars.push_back(1);
                detailed_stars.push_back(10);
                detailed_stars.push_back(18);
                detailed_stars.push_back(22);
    cout << "new alpha new luck" << endl;
    /*
                for(int l=0;l<detailed_stars.size();l++)
                {
                    string detailed_star = to_string(detailed_stars[l]);
                    output_path = "./Plotting/Results/Profiles/star";
                    output_path.append(to_string(l).c_str());
                    output_path.append("/TOV_output_");
                    output_path.append(theory_ID);

                    p_rho_profiles.push_back(fopen(output_path.c_str(),"w"));
                    rewind(p_rho_profiles[l]);
                }
    */
                while(rho_center_now < rho_center_max)
                {
                    ccount++;
                    pair<double, double> Mass_Radius = tov_integrate(rho_center_now);

                    Mass = Mass_Radius.first;
                    Radius = Mass_Radius.second;

                    Masses.push_back(Mass);
                    Radi.push_back(Radius);

                    cout << "Mass at Star# " << ccount << " = "<< Mass << endl;
                    cout << "rho at star# " << ccount << " = " << rho_center_now << endl;
                    cout << "Radius at star# " << ccount << " = " << Radius << endl << endl;

                    //Write Result on File "TOV_output"
                    if(alpha_iterate==0 || beta_iterate==0)
                    {
                        fprintf(ifp,"%.10lf %.10lf %.10lf\n",rho_center_now,Mass,Radius);
                        fprintf(ifp2,"%.10lf %.10lf %.10lf\n",rho_center_now,Mass,Radius);
                    }

                    rho_center_now = rho_center_now*1.2;
                    if(Mass>Maxmass)
                    {
                        Maxmass = Mass;
                        Maxmass_Radius = Radius;
                    }
                    cout << "rho_center_now= " << rho_center_now << endl;

                }
                rewind(ifp);
                rewind(ifp2);

                ///////
                // Save Maxmass of this alpha
                if(alpha_iterate==1 && beta_iterate==0){ fprintf(ifp_maxmass_alpha,"%.10lf %.10lf %.10lf\n",alpha,Maxmass,Maxmass_Radius); }

                ///////
                // Save Maxmass of this alpha
                if(alpha_iterate==0 && beta_iterate==1){ fprintf(ifp_maxmass_beta,"%.10lf %.10lf %.10lf\n",beta,Maxmass,Maxmass_Radius); }

                ///////
                // Save Maxmass of this alpha
                if(alpha_iterate==1 && beta_iterate==1){ fprintf(ifp_maxmass_alpha_beta,"%.10lf %.10lf %.10lf %.10lf\n",alpha,beta,Maxmass,Maxmass_Radius);}
                }



            /////////////////////////
            //PLOT TOV-RESULTS: R,M
            ////////////////////////
            //plotting_function(File ifp, string file_name, string title, string x_label, string y_label, bool x_log, bool y_log, int x_eq_col_num, int y_eq_col_num)
            if(alpha_iterate==0 && beta_iterate==0){ plotting_function("TOV_output", "R to M", "R[cm]", "M/M_O");}
            fclose(ifp);
            fclose(ifp2);
            cout << "Maxmass= " << Maxmass << endl;
            Maxmass = 0;
            Maxmass_Radius = 0;

            ccount = 0;
            beta=beta+beta_step;
            Rp = 1.0/alpha;
            Rq = 1.0/beta;
            if(beta_iterate==0 || fRQ_theory!=1){break;}
            cout << "===================================" << endl;

        } // Beta-close

        ccount = 0;
        alpha=alpha+alpha_step;
        Rp = 1.0/alpha;
        Rq = 1.0/beta;
        if(alpha_iterate==0){break;}
        cout << "===================================" << endl;

    } // ALPHA-close
    rewind(ifp_maxmass_alpha);
    rewind(ifp_maxmass_beta);
    rewind(ifp_maxmass_alpha_beta);
}

/////////////////////END MAIN

