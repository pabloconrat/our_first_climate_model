#ifndef __reptran_thermal_h
#define __reptran_thermal_h

#if defined (__cplusplus)
extern "C" {
#endif

void read_tau (const char *reducedLkpPath, int nLev, double *p, double *T, double *H20_VMR,
	       double *CO2_VMR, double *O3_VMR, double *N2O_VMR, double *CO_VMR, double *CH4_VMR,
	       double *O2_VMR, double *HNO3_VMR, double *N2_VMR,
	       double ***tau, double **wvl, double **weight, int *nWvl);

#if defined (__cplusplus)
}
#endif

#endif
