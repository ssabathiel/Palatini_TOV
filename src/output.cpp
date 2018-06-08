#include "header/eos_filetoarray.h"

using namespace std;

////////////////////////
// WRITES RESULTS ON FILES
////////////////////////

//Get output (collumns not seperated by "," as in CSV, but by space " ") to array
void output_to_arrays(const char* i_file_name,  vector<double> &rhos, vector<double> &masses, vector<double> &radii )
{
        // removed all the windowing stuff.
    ifstream in;
    in.open(i_file_name);

    double tmp_rho;
    double tmp_mass;
    double tmp_radius;

    if (in.is_open())
    {
        // read each column into an array
        for (int i = 0; i < 364; i++)
        { // what do you do if the file ends early or contains bad information?
            in >> tmp_rho >>
                  tmp_mass >>
                  tmp_radius;
            rhos.push_back(tmp_rho);
            masses.push_back(tmp_mass);
            radii.push_back(tmp_radius);
        }
        in.close();
    }

    else
    {
        cout << "error reading file" << endl;
        return;
    }
}

string DoublePowToString(double Rpp)
{
    // TODO:
    // Rephrase this part into: Base_val_shift, Base_pos_shift, Pot_val_shift, Pot_pos_shift. And go serial +/- by checking each condition
    cout << "Rpp== " << Rpp << endl;

    int pot_val_shift=0;
    int pot_pos_shift=0;

    string coeff_sign;
    int shift=0;
    if(Rpp<0)
    {
        coeff_sign="";
        shift = 1;

    }
    else
    {
        coeff_sign="+";
    }

    string pow_sign;
    int shift2=0;
    int shift3=0;
    int shift5=0;
    if(log10(abs(Rpp))<0)
    {
        pow_sign="";
        shift2=0;
        pot_pos_shift += 1;
        if(abs(log10(abs(Rpp)))>9)
        {
            pot_pos_shift += 1;
            pot_val_shift += 0;

        }
        else
        {
            pot_val_shift += 0;
            pot_pos_shift -= 1;
        }

    }
    else
    {
        pow_sign="+";
        if(abs(log10(abs(Rpp)))>9.999)
        {
            pot_pos_shift += 1;
        }
        else
        {
            pot_pos_shift -= 1;
        }
        shift3=1;
        shift2=0;
    }

    if(abs(log10(abs(Rpp)))>9)
    {
        cout << "Yes  bigger 9" << endl;
        shift5=1;
    }

    int shift4=0;
    if(abs(log10(abs(Rpp)))>11-shift3)
    {
        shift4=1;
    }

    cout << "log10(abs(Rpp))= " << log10(abs(Rpp)) << endl;
    char Rp_value[4];
    char Rp_pot[3+pot_pos_shift];//shift3+shift4+shift5];

    memset(Rp_value, 0, sizeof(Rp_value));
    snprintf(Rp_value, sizeof(Rp_value), "%g", Rpp);
    std::string strObj4(Rp_value);

    memset(Rp_pot, 0, sizeof(Rp_pot));
    snprintf(Rp_pot, sizeof(Rp_pot), "%g", log10(abs(Rpp)) - pot_val_shift);//-shift2+shift4+shift5 );
    std::string strObj2(Rp_pot);




    cout << Rp_value[0] << endl;
    cout << Rp_value[1] << endl;
    cout << Rp_value[2] << endl;

    if(Rp_value[1+shift]=='e')
    {
        Rp_value[1+shift]= '.';
        Rp_value[2+shift] = '0';
    }

    string Rp_string=coeff_sign;
    Rp_string.append(Rp_value);
    Rp_string.append("*10^");
    Rp_string.append(pow_sign);
    Rp_string.append(Rp_pot);


    return Rp_string;

}

string create_theory_ID_wo_alpha_beta()
{
    string theory_ID;
    string EOS_ID_str;
    string Rp_string;
    string Rq_string;


    if(analytical_EOS==1){EOS_ID_str="_SLY";}
    else if(tabular_EOS==1){EOS_ID_str="_TAB";}
    else if(ply_EOS==1){EOS_ID_str="_PLY";}
    else if(fps_EOS==1){EOS_ID_str="_FPS";}
    else if(ap4_EOS==1){EOS_ID_str="_AP4";}

    if(th==0)
    {
        theory_ID  = "GR";
        theory_ID.append(EOS_ID_str);
    }
    if(th==1 || th==4)
    {
        theory_ID  = "fR";
        theory_ID.append(EOS_ID_str);

    }
    if(th==2)
    {
        theory_ID  = "fRQ";
        theory_ID.append(EOS_ID_str);
    }

    return theory_ID;
}



string create_theory_ID()
{
    string theory_ID;
    string EOS_ID_str;
    string Rp_string;
    string Rq_string;

    Rp_string = DoublePowToString(alpha);
    Rq_string = DoublePowToString(beta);

    if(analytical_EOS==1){EOS_ID_str="_SLY";}
    else if(tabular_EOS==1){EOS_ID_str="_TAB";}
    else if(ply_EOS==1){EOS_ID_str="_PLY";}
    else if(fps_EOS==1){EOS_ID_str="_FPS";}
    else if(ap4_EOS==1){EOS_ID_str="_AP4";}

    if(th==0)
    {
        theory_ID  = "GR";
        theory_ID.append(EOS_ID_str);
    }
    if(th==1 || th==4)
    {
        theory_ID  = "fR_Rp_";
        theory_ID.append(Rp_string);
        theory_ID.append(EOS_ID_str);

    }
    if(th==2)
    {
        theory_ID  = "fRQ_Rp_";
        theory_ID.append(Rp_string);
        theory_ID.append("_Rq_");
        theory_ID.append(Rq_string);
        theory_ID.append(EOS_ID_str);
    }

    return theory_ID;
}




