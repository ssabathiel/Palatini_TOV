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
bool fRQ_theory;

bool analytical_EOS;
bool ply_EOS;
bool tabular_EOS;

double Rp;
double Rq;
int th;
int num_an;

vector<rho_of_p_functions> rho_of_p;
vector<p_of_rho_functions> p_of_rho;
vector<drho_dp_functions> drho_dp;
vector<ddrho_dpp_functions> ddrho_dPP;

int main()
{

    configure();
    rho_of_p.push_back(rho_of_p_num);
    rho_of_p.push_back(rho_of_p_analytical);
    rho_of_p.push_back(rho_of_p_analytical_ply);

    p_of_rho.push_back(p_of_rho_num);
    p_of_rho.push_back(p_of_rho_analytical);
    p_of_rho.push_back(p_of_rho_analytical_ply);

    drho_dp.push_back(drho_dp_num);
    drho_dp.push_back(drho_dp_analytical);
    drho_dp.push_back(drho_dp_analytical_ply);

    ddrho_dPP.push_back(ddrho_dPP_num);
    ddrho_dPP.push_back(ddrho_dPP_analytical);
    ddrho_dPP.push_back(ddrho_dPP_analytical_ply);


    if(GR_theory==1){th=0;}
    else if(fR_theory==1){th=1;}
    else if(fRQ_theory==1){th=2;}
    else if(fR_lim_theory==1){th=3;}

    if(analytical_EOS==1){num_an=0;}
    else if(tabular_EOS==1){num_an=1;}
    else if(ply_EOS==1){num_an=2;}

    char* test_dat = "EOS_data/Sly_csv";    //test_data.dat
    create_spline(test_dat);

    ///////////////////////////////////////
    //LOOP OVER RHO_CENTERS if multiple stars, other wise INTEGRATE ONLY ONCE
    ///////////////////////////////////////
    double Mass;
    double Radius;
    vector<double> Masses;
    vector<double> Radi;

    double rho_center_now = rho_center;
    double rho_center_max = rho_center_max2;

    string theory_ID = create_theory_ID();

    string output_path = "./Plotting/Results/TOV_output_";
    output_path.append(theory_ID);
    FILE *ifp=fopen(output_path.c_str(),"w");
    rewind(ifp);
    FILE *ifp2=fopen("TOV_output","w");
    rewind(ifp2);


    vector<int> detailed_stars;
    detailed_stars.push_back(7);
    detailed_stars.push_back(41);
    detailed_stars.push_back(48);
    detailed_stars.push_back(49);


    for(int l=0;l<detailed_stars.size();l++)
    {
        string detailed_star = to_string(detailed_stars[l]);
        output_path = "./Plotting/Results/Profiles/p_of_rho_ccount_";
        output_path.append(detailed_star);
        p_rho_profiles.push_back(fopen(output_path.c_str(),"w"));
        rewind(p_rho_profiles[l]);
    }


    while(rho_center_now < rho_center_max){
        ccount++;
        pair<double, double> Mass_Radius = tov_integrate(rho_center_now);
        Mass = Mass_Radius.first;
        Radius = Mass_Radius.second;
        if(Mass < 3000.0)          //Due to numerical instabilities
        {
        Masses.push_back(Mass);
        Radi.push_back(Radius);
        }
        cout << "Mass at Star# " << ccount << " = "<< Mass << endl;
        cout << "rho at star# " << ccount << " = " << rho_center_now << endl;
        cout << "Radius at star# " << ccount << " = " << Radius << endl << endl;

        //Write Result on File "TOV_output"
        //if(Mass < 3.0 && Mass >0)          //Due to numerical instabilities
        if(Mass >0)          //Due to numerical instabilities
        {
            fprintf(ifp,"%.10lf %.10lf %.10lf\n",rho_center_now,Mass,Radius);
            fprintf(ifp2,"%.10lf %.10lf %.10lf\n",rho_center_now,Mass,Radius);
        }

        rho_center_now = rho_center_now*1.2;
        cout << "rho_center_now= " << rho_center_now << endl;

    }
    rewind(ifp);
    rewind(ifp2);



    /////////////////////////
    //PLOT TOV-RESULTS: R,M
    ////////////////////////
    //plotting_function(File ifp, string file_name, string title, string x_label, string y_label, bool x_log, bool y_log, int x_eq_col_num, int y_eq_col_num)
    plotting_function("TOV_output", "R to M", "R[cm]", "M/M_O");
    fclose(ifp);
    fclose(ifp2);


}

/////////////////////END MAIN

