#include <stdio.h>
#include "ascii.h"

void do_nothing (double *wvl, double **tau, int nwvl, int nlyr)
{
  // do nothing, only show how 1D and 2D arrays are passed to a function

  return;
}



int main (int argc, char **argv)
{
  int nwvl=0, nlyr=0;
  int status=0;
    
  char tauCO2filename[FILENAME_MAX]="./lbl.arts/lbl.co2.asc";

  double *wvl=NULL;        // 1D array
  double **tauCO2=NULL;    // 2D array 

  /* read CO2 optical thickness profile */
  status = ASCII_file2xy2D (tauCO2filename,   
			    &nwvl, &nlyr, 
			    &wvl, &tauCO2);
  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, tauCO2filename);
    return status;
  }

  fprintf (stderr, " ... read %d wavelengths and %d layers from %s\n",
	   nwvl, nlyr, tauCO2filename);

  do_nothing (wvl, tauCO2, nwvl, nlyr);
  
  return 0;  
}
