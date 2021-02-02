#ifndef __reptran_thermal_h
#define __reptran_thermal_h

void read_tau (const char *reducedLkpPath, int nLev, std::vector<double> &plevel, std::vector<double> &Tvector, double *H20_VMR,
	       double *CO2_VMR, double *O3_VMR, double *N2O_VMR, double *CO_VMR, double *CH4_VMR,
	       double *O2_VMR, double *HNO3_VMR, double *N2_VMR,
	       double ***tau, double **wvl, double **weight, int *nWvl,
	       int T_at_Lev);

#endif
