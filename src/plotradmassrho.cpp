#include "header/plotradmassrho.h"
#include <fstream>
using namespace std;

//////////////////////////////
//PLOTTTING FUNCTION FROM FILE
//////////////////////////////


void plotting_function(char * file_name, char * title, char * x_label, char * y_label)
{

    //char path[]="./";
    char path[]="/home/silvester/fRQ_Projectt/build-fRQ_TOV2-Desktop-Debug/";

    FILE *gnuplot = popen("/usr/bin/gnuplot -persist","w");
    string fily = "TOV_output";

    fprintf(gnuplot,"set multiplot\n");
    fprintf(gnuplot,"set origin 0.0,0.0\n");
    fprintf(gnuplot,"set size 0.98,0.48\n");
    fprintf(gnuplot,"set logscale x\n");
    fprintf(gnuplot,"set ylabel \"M/M_O\"\n");
    fprintf(gnuplot,"set xlabel \"{/Symbol r} (g/cm^{-3})\"\n");
    fprintf(gnuplot,"set style data linespoints\n");

    fprintf(gnuplot,"plot \"%s%s\" using 1:2 notitle lt 3\n",path,fily.c_str());
    fprintf(gnuplot,"set origin 0.00,0.50\n");
    fprintf(gnuplot,"set size 0.98,0.48\n");
    fprintf(gnuplot,"set ylabel \"%s\"\n",y_label);
    fprintf(gnuplot,"set xlabel \"%s\"\n",x_label);
    //fprintf(gnuplot,"plot \"%sTOV_output_fR\" using 3:2 notitle lt 3\n",path);
    fprintf(gnuplot,"plot \"%s%s\" using 3:2 notitle lt 3\n",path,fily.c_str());

    fprintf(gnuplot,"exit \n");
    pclose(gnuplot);


}


