#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "ascii.h"
#include "repwvl_thermal.h"

#define SIGMA 5.67E-8
#define G 9.8065
#define R 287.0
#define CP 1004.0

#define K_BOLTZMANN 1.3806503E-23
#define H_PLANCK 6.626068E-34
#define C_LIGHT 299792458

double planck (double T, double wavelength)
{
  wavelength*=1e-9;  /* convert from nm to m */

  return 2.0*H_PLANCK*C_LIGHT*C_LIGHT/wavelength/wavelength/wavelength/wavelength/wavelength/(exp(H_PLANCK*C_LIGHT/wavelength/K_BOLTZMANN/T)-1.0)/1.0e9;  /* W/(m2 nm sterad) */
}


int main (int argc, char **argv)
{
  int nlev=0, nwvl=0, status = 0;


  double *plev=NULL, *plevPa=NULL, *Tlev=NULL, *rholev=NULL;
  double *wvl=NULL;
  double *tmp=NULL;
  double Bg=0;

  double *H2O_VMR=NULL;
  double *CO2_VMR=NULL;
  double *O3_VMR=NULL;
  double *N2O_VMR=NULL;
  double *CO_VMR=NULL;
  double *CH4_VMR=NULL;
  double *O2_VMR=NULL;
  double *HNO3_VMR=NULL;
  double *N2_VMR=NULL;
  
  /* first read atmospheric profile */
  status = read_9c_file ("test.atm", &tmp, &plev, &Tlev, &rholev, &H2O_VMR, &O3_VMR, &CO2_VMR, &CH4_VMR, &N2O_VMR, &nlev);
  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, "test.atm");
    return status;
  }
    
  O2_VMR   = calloc (nlev, sizeof(double));
  N2_VMR   = calloc (nlev, sizeof(double));
  HNO3_VMR = calloc (nlev, sizeof(double));
  CO_VMR   = calloc (nlev, sizeof(double));

  plevPa   = calloc (nlev, sizeof(double));
  
  for (int ilev=0; ilev<nlev; ilev++) {
    /* convert from ppm to absolute concentrations */
    H2O_VMR [ilev] *= 1E-6;
    O3_VMR  [ilev] *= 1E-6;
    CO2_VMR [ilev] *= 1E-6;
    CH4_VMR [ilev] *= 1E-6;
    N2O_VMR [ilev] *= 1E-6;
    
    /* define missing species */
    O2_VMR  [ilev] = 0.2095;
    N2_VMR  [ilev] = 0.7808;
    HNO3_VMR[ilev] = 0.0;
    CO_VMR  [ilev] = 0.0;

    /* convert pressure from hPa to Pa */
    plevPa  [ilev] = plev[ilev] * 100.0;
  }

  /* obtain spectral absorption coefficients */
  
  double **tau = NULL;
  double *weight = NULL;

  read_tau ("./Reduced100Forcing.nc", nlev, plevPa, Tlev,
	    H2O_VMR, CO2_VMR, O3_VMR, N2O_VMR, CO_VMR, CH4_VMR, O2_VMR, HNO3_VMR, N2_VMR,
	    &tau, &wvl, &weight, &nwvl, 1);
  
  for (int iv=0; iv<nwvl; iv++)
    Bg += weight[iv] * planck (Tlev[nlev-1], wvl[iv]);
  
  fprintf (stderr, "Bg = %f\n", M_PI*Bg);

  return 0;  
}



