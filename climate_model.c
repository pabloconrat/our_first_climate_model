//Climate Model Version 5:21 PM 30/01/2017


 
//Import librarues
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ascii.h"
#include <unistd.h>
#include "gnuplot_i.h"
#include "fpda_rrtm_lw.h"
#include "fpda_rrtm_sw.h"
  
//global variables
#define c_p 1003.
#define p_ground 100000.
#define wavelength_short 100	//shortest wavelength for discrete wavlength calculation
#define wavelength_long 103	//longest wavelength for descrete wavelength calculation
#define delta_wavelength 1	//descrete wavelength between wavelength_short and wavelength_long
#define delta_angle 45		//0 degree, 90 degree and 8 angles in between
#define g 9.81
#define M_PI 3.14159265359
  
 
//--------------------------------- prototype ------------------------------

double StephanBoltzmann (double Temperature);

int To_T_pot (double *T, double *T_pot, int layer, double *p);

int T_pot_To_T (double *T, double *T_pot, int layer, double *p);

int WriteToFile (char *WordText);

int Sort_Layers (double *T_pot, int layer);

int cmpfunc (const void *x, const void *y);

int HydrostaticEquation (int layer, double *z, double *T, double *p);

double integrate_B_Planck (double T, double lambda_min, double lambda_max);

int Adiabatic_T_Profile (int layer, double T_ground, double *T, double *p);

int Calculate_layer_p (double *p, int layer);

 
int Thermal_RT (double *T, int layer, double *E_up, double *E_down,
		   int nbands, double *band_lbound, double *band_ubound,
		   double **wgt_lw, double **dtau_mol_lw,
		   double *omega0_cloud_lw, double *tau_cloud);

int Schwarzschild_RT (double *B, double *tau, int layer, double Bg,
		       double *E_down, double *E_up);

 
void eddington_v2 (double dtau, double g_asymm, double omega0, double mu0,
		      double *t, double *r, double *rdir, double *sdir,
		      double *tdir);

int Visual_RT (int layer, double *E_up_vis, double *E_down_vis, double *Edir,
		int nbands, double **dtau_mol_sw, double **dtau_ray_sw,
		double *wgt_sw, double mu0, double albedo,
		double *omega0_cloud_sw, double *tau_cloud, double *g_cloud);

void TwoStream_RT (int layer, double *tau, double *w0, double *g_asymm,
		    double S0, double mu, double *E_up_spec,
		    double *E_down_spec, double *Edir_spec, double albedo);

 
int Define_Cloud (int layer, double tau_cloud_init, double *tau_cloud,
		     int cloud_layer);

int DeltaScaling (int il, int ib, double *omega0_cloud, double *g_cloud,
		   double *tau_cloud, double omega0_cloud_star,
		   double g_cloud_star, double tau_cloud_star);

 
int Rel_hum_To_vmr (double *vmr, double *esat, double *p, double *rel_hum,
		       int layer);

int Magnus_formula (double *esat, double *T, int layer);

int Rel_hum_To_vmr (double *vmr, double *esat, double *p, double *rel_hum,
		     int layer);

//--------------------------------- prototype End------------------------------
  
 
 
//--------------------------------------------- Functions ------------------------------------------------------
  
//StephanBoltzmann Temperature to Energy
  double
StephanBoltzmann (double Temperature)
{
  
double sigma = 5.67e-8;
  
return sigma * pow (Temperature, 4);

}


 
 
 
//Calulate the potential Temperatur from Temperatur and height
  int
To_T_pot (double *T, double *T_pot, int layer, double *p)
{
  
double R = 287.;
  
for (int i = 0; i < layer; i++)
    {
      
T_pot[i] = T[i] * pow ((p_ground / (p[i])), (R / c_p));
    
} 
 
return 0;

}


 
 
 
//Calulate the Temperatur from potential Temperatur and height
  int
T_pot_To_T (double *T, double *T_pot, int layer, double *p)
{
  
double R = 287.;
  
for (int i = 0; i < layer; i++)
    {
      
T[i] = T_pot[i] * pow ((p_ground / (p[i])), (-R / c_p));
    
} 
 
return 0;

}


 
 
 
//Compare function for qsort function
  int
cmpfunc (const void *x, const void *y)
{
  
double xx = *(double *) x, yy = *(double *) y;
  
if (xx < yy)
    return 1;
  
if (xx > yy)
    return -1;
  
return 0;

}


 
 
 
//Advection functions, sorting the layers by potential Temperature
  int
Sort_Layers (double *T_pot, int layer)
{
  
qsort (T_pot, layer, sizeof (double), cmpfunc);
  
return 0;

}


 
 
 
//Calculates the height in meters from the pressure
  int
HydrostaticEquation (int layer, double *z, double *T, double *p)
{
  
double R = 287.;
  
double z_culevel = 0.;
  
double delta_z = 0.;
  
double delta_p = p_ground / layer;
  
for (int i = layer - 1; i >= 0; i--)
    {
      
delta_z = R * T[i] * delta_p / (g * (p[i]));
      
z[i] = z_culevel + 0.5 * delta_z;
      
z_culevel += delta_z;
    
} 
return 0;

}


 
 
double
integrate_B_Planck (double T, double lambda_min, double lambda_max)
{
  
double precision;
  
if (lambda_min < 20e-6)
    {
      
precision = 100e-9;
    
}
  else
    {
      
precision = 500e-9;
    
}
  
int nplanck = 1 + ((lambda_max - lambda_min) / precision);	//---------------------------------- integration steps wavelength ----------------------------------
  double B_planck = 0.;		//initialize
  double c_light = 299792458.;	//m/s
  double h_planck = 6.626070e-34;	//J*s
  double k_B = 1.38064852e-23;	//J/K   
  
for (int i = 0; i < nplanck; i++)
    {
      
double wavelength =
	lambda_min + ((i + 0.5) * ((lambda_max - lambda_min) / nplanck));
      
double B_Planck_singleband =
	((2 * h_planck * pow (c_light, 2)) /
	 (pow (wavelength, 5) *
	  (exp ((h_planck * c_light) / (k_B * wavelength * T)) - 1)));
      
 
B_planck +=
	((lambda_max - lambda_min) / nplanck) * B_Planck_singleband;
    
 
} 
return B_planck;

}


 
 
 
int
Adiabatic_T_Profile (int layer, double T_ground, double *T, double *p)
{
  
    //Define adiabatic temperature profile, Initial conditions
  double R = 287.;
  
double tempgrad = 0.0065;	// K/m
  double delta_p = p_ground / layer;
  
double T_lower = T_ground;
  
 
double T_upper =
    T_lower / (1 + (tempgrad * R * delta_p) / (g * p_ground));
  
T[layer - 1] = 0.5 * (T_lower + T_upper);
  
for (int i = layer - 1; i > 0; i--)
    {
      
T_lower = T[i] / (1 + (tempgrad * R * delta_p) / (2 * g * p[i]));
      
T_upper = T_lower / (1 + (tempgrad * R * delta_p) / (g * p[i]));
      
T[i - 1] = 0.5 * (T_lower + T_upper);
    
} 
return 0;

}


 
 
 
//calculate discrete layer pressures
  int
Calculate_layer_p (double *p, int layer)
{
  
for (int i = 0; i < layer; i++)
    {
      
p[i] = i * (p_ground / layer) + 0.5 * (p_ground / layer);
    
} 
return 0;

}


 
 
//bla
//Calculates the radiative transfer for given B and tau; one wavelength
  int
Schwarzschild_RT (double *B, double *tau, int layer, double Bg,
		  double *E_down, double *E_up)
{
  
    //variables
  double L_down[layer + 1];	//definiert als level
  double L_up[layer + 1];
  
double tau_total[layer];
  
 
 
    //Reset von E and calculate tau_total for all layers
    for (int i = 0; i < layer + 1; i++)
    {
      
 
 
E_up[i] = 0;
      
E_down[i] = 0;
    
} 
 
 
 
 
    //loop over descrete angles
  int nangles = 8;		//---------------------------------- integration steps angles ----------------------------------
  for (int iangle = 0; iangle < nangles; iangle++)
    {
      
double mu = (1. / (2. * nangles)) + ((float) iangle / (float) nangles);
      
 
 
	//------------ Calculate from top to bottum get E_down ---------------
	
L_down[0] = 0;		//initialize
      
for (int i = 0; i < layer; i++)
	{
	  
	    //transmission plus emmision
	    L_down[i + 1] =
	    L_down[i] * exp (-tau[i] / mu) + B[i] * (1 - exp (-tau[i] / mu));
      
} 
 
 
for (int i = 0; i < layer + 1; i++)
	{
	  
E_down[i] += L_down[i] * mu * 2 * M_PI * (1. / nangles);
	
 
} 
 
	//------------ Calculate from top to bottom get E_down END---------------
	
 
 
	//------------ Calculate from top to bottom get E_up ---------------
	
L_up[layer] = Bg;	//initialize Background radiation from ground
      
for (int i = layer; i > 0; i--)
	{
	  
L_up[i - 1] =
	    L_up[i] * exp (-tau[i - 1] / mu) + B[i - 1] * (1 -
							   exp (-tau[i - 1] /
								mu));
	  
	    //                  transmission           +        emmision
      } 
 
for (int i = layer; i >= 0; i--)
	{
	  
E_up[i] += L_up[i] * mu * 2 * M_PI * (1. / nangles);
	
} 
	//------------ Calculate from bottom to top get E_down up---------------
    
 
}			//end of angle loop
  
return 0;

}


 
 
 
 
 
 
//Radiative read in 
  int
Thermal_RT (double *T, int layer, double *E_up, double *E_down, int nbands,
	    double *band_lbound, double *band_ubound, double **wgt_lw,
	    double **dtau_mol_lw, double *omega0_cloud_lw, double *tau_cloud)
{
  
 
double E_up_spec[layer + 1];	// auf level definiert
  double E_down_spec[layer + 1];
  
double B[layer], Bg;
  
 
 
    //reset E_up and E_down
    for (int i = 0; i < layer + 1; i++)
    {
      
E_up[i] = 0;
      
E_down[i] = 0;
    
} 
 
    //loop over the given wavelength array
    for (int ib = 0; ib < nbands; ib++)
    {
      
 
	//inital conditions for each wavelength
      double wavelength_start = 1e-2 / band_ubound[ib];
      
double wavelength_stop = 1e-2 / band_lbound[ib];
      
Bg =
	integrate_B_Planck (T[layer - 1], wavelength_start,
			    wavelength_stop) * wgt_lw[ib][layer - 1];
      
 
double tau_lw[layer];
      
for (int il = 0; il < layer; il++)
	{
	  
 
	    //Calculate B  for each layer with a given wavelength spec
	    B[il] =
	    integrate_B_Planck (T[il], wavelength_start,
				wavelength_stop) * wgt_lw[ib][il];
	  
 
	    //Calculate tau for all layers at given band
	    tau_lw[il] = dtau_mol_lw[ib][il] + (1 - omega0_cloud_lw[ib]) * tau_cloud[il];	//only absorption of the cloud ==> (1-omage0)
	
} 
 
 
	// thermal radiative transfer               
	//Radiative transfer over angles with given B[layer] and tau_lw[layer]
	Schwarzschild_RT (B, tau_lw, layer, Bg, E_down_spec, E_up_spec);
      
 
	//Sum up all E over all wavelengths
	for (int il = 0; il < layer + 1; il++)
	{
	  
E_up[il] += E_up_spec[il];
	  
E_down[il] += E_down_spec[il];
    
} 
 
 
}			//End of wavelength loop     
  
return 0;

}


 
 
int
Visual_RT (int layer, double *E_up_vis, double *E_down_vis, double *Edir,
	   int nbands, double **dtau_mol_sw, double **dtau_ray_sw,
	   double *wgt_sw, double mu0, double albedo, double *omega0_cloud_sw,
	   double *tau_cloud, double *g_cloud)
{
  
 
double E_up_spec[layer + 1];	// auf level definiert
  double E_down_spec[layer + 1];
  
double Edir_spec[layer + 1];
  
 
//for(int ib = 0; ib<100; ib++){
//              for(int il = 0; il<layer; il++){
    
 
    //reset E_up and E_down
    for (int i = 0; i < layer + 1; i++)
    {
      
E_up_vis[i] = 0;
      
E_down_vis[i] = 0;
      
Edir[i] = 0;
    
} 
 
 
    //loop over the given wavelength array
    for (int ib = 0; ib < nbands; ib++)
    {
      
 
 
double tau_sw[layer];
      
double w0[layer];
      
double g_asymm[layer];
      
 
for (int il = 0; il < layer; il++)
	{
	  
 
double f = g_cloud[ib] * g_cloud[ib];
	  
double tau_cloud_star =
	    (1 - f * omega0_cloud_sw[ib]) * tau_cloud[il];
	  
double g_cloud_star = (g_cloud[ib] - f) / (1 - f);
	  
double omega0_cloud_star =
	    ((1 - f) * omega0_cloud_sw[ib]) / (1 - f * omega0_cloud_sw[ib]);
	  
 
tau_sw[il] = dtau_mol_sw[ib][il] + dtau_ray_sw[ib][il] + tau_cloud_star;	//hier cloud scat und abs
	  w0[il] = (dtau_ray_sw[ib][il] + omega0_cloud_star * tau_cloud_star) / tau_sw[il];	//single scatter albedo
	  g_asymm[il] = (0 * dtau_ray_sw[ib][il] + omega0_cloud_star * tau_cloud_star * g_cloud_star) / (w0[il] * tau_sw[il]);	//forward to backscatter ratio; 0 means same amount forward and backward ; changes with clouds
	
} 
 
	//radiative transfer in visual range
	TwoStream_RT (layer, tau_sw, w0, g_asymm, wgt_sw[ib], mu0, E_up_spec,
		      E_down_spec, Edir_spec, albedo);
      
 
 
	//Sum up all E over all wavelengths
	for (int il = 0; il < layer + 1; il++)
	{
	  
E_up_vis[il] += E_up_spec[il];
	  
E_down_vis[il] += E_down_spec[il];
	  
Edir[il] += Edir_spec[il];
    
} 
 
 
}			//End of wavelength loop     
  
 
return 0;

}


 
 
 
 
void
TwoStream_RT (int layer, double *tau, double *w0, double *g_asymm, double S0,
	      double mu, double *E_up_spec, double *E_down_spec, double *Edir,
	      double albedo)
{
  
 
double t[layer];
  
double r[layer];
  
double rdir[layer];
  
double sdir[layer];
  
double tdir[layer];
  
double Tdir[layer];
  
double Rdir[layer];
  
double Sdir[layer];
  
double T[layer];
  
double R[layer];
  
 
 
for (int il = 0; il < layer; il++)
    {
      
eddington_v2 (tau[il], g_asymm[il], w0[il], mu, &(t[il]), &(r[il]),
		     &(rdir[il]), &(sdir[il]), &(tdir[il]));
    
} 
 
    //Randbedingungen
    Edir[0] = S0 * mu;		//ToA
  E_down_spec[0] = 0;		//reset E_down_spec
  E_up_spec[0] = 0;
  
Tdir[0] = tdir[0];
  
Sdir[0] = sdir[0];
  
R[0] = r[0];
  
T[0] = t[0];
  
 
 
    //--------------Doubling Adding Scheme------------------
    
    //Calculate R, T, Tdir, Sdir, Edir
    for (int il = 0; il < layer - 1; il++)
    {
      
R[il + 1] =
	r[il + 1] +
	((R[il] * t[il + 1] * t[il + 1]) / (1 - R[il] * r[il + 1]));
      
T[il + 1] = (T[il] * t[il + 1]) / (1 - R[il] * r[il + 1]);
      
	//printf("il=%i g=%f tau=%f w0=%f r=%f\n",il, g_asymm[il],tau[il],w0[il],r[il]);        
	//Calculate T_dir with a product of of all t_dir from layer 0 to layer
	Tdir[il + 1] = Tdir[il] * tdir[il + 1];
      
	//Calculate Sdir from ...? Formal Script 13 doubling adding
	Sdir[il + 1] =
	((t[il + 1] * Sdir[il] +
	  Tdir[il] * rdir[il + 1] * R[il] * t[il + 1]) / (1 - R[il] * r[il +
									1])) +
	(Tdir[il] * sdir[il + 1]);
    
 
} 
 
 
    //Calculate Edir from lambert beer 
    for (int il = 0; il < layer; il++)
    {
      
Edir[il + 1] = Edir[il] * tdir[il];
    
} 
 
 
 
 
 
 
 
    //Calculate irradience at the surface
    E_down_spec[layer] =
    ((Sdir[layer - 1] + Tdir[layer - 1] * R[layer - 1] * albedo) / (1 -
								    R[layer -
								      1] *
								    albedo)) *
    Edir[0];
  
E_up_spec[layer] =
    albedo * (E_down_spec[layer] + Tdir[layer - 1] * Edir[0]);
  
 
 
 
 
 
//Calculate irradience for each layer starting from the bottom
    for (int il = layer; il > 1; il--)
    {
      
	//                                                  no source term script page 7                                                                    source term script page 11
	E_down_spec[il - 1] =
	((T[il - 2] * E_down_spec[0] +
	  R[il - 2] * t[il - 1] * E_up_spec[il]) / (1 - R[il - 2] * r[il -
								      1])) +
	((Edir[0] * Sdir[il - 2] +
	  Edir[il - 1] * rdir[il - 1] * R[il - 2]) / (1 - R[il - 2] * r[il -
									1]));
      
 
E_up_spec[il - 1] =
	((t[il - 1] * E_up_spec[il] +
	  T[il - 2] * r[il - 1] * E_down_spec[0]) / (1 - R[il - 2] * r[il -
								       1])) +
	((Edir[0] * Sdir[il - 2] * r[il - 1] +
	  Edir[il - 1] * rdir[il - 1]) / (1 - R[il - 2] * r[il - 1]));
      
//printf("%f\n",Rdir[il-1]);
    } 
 
 
    //Is this correct???? ->
    
    //diffuse and direct componet at top of atmosphere
    E_up_spec[0] = t[0] * E_up_spec[1] + rdir[0] * Edir[0];

 
 
} 
 
 
 
void

eddington_v2 (double dtau, double g_asymm, double omega0, double mu0,
	      double *t, double *r, double *rdir, double *sdir, double *tdir)
{
  
 
    /* calculate layer properties t, r, rdir, sdir, and tdir from            */ 
    /* layer optical thickness dtau, asymmetry parameter g,                  */ 
    /* single scattering albedo omega0, and cosine of solar zenith angle mu0 */ 
  double alpha1 = 0, alpha2 = 0, alpha3 = 0, alpha4 = 0, alpha5 = 0, alpha6 =
    0;
  
double a11 = 0, a12 = 0, a13 = 0, a23 = 0, a33 = 0;
  
double lambda = 0, b = 0, A = 0;
  
double denom = 0;
  
 
    /* first, avoid omega0=1 because of instability */ 
    if (omega0 > 0.999999)
    
omega0 = 0.999999;
  
 
if (dtau > 100)
    
dtau = 100;
  
 
alpha1 = (1.0 - omega0) + 0.75 * (1.0 - omega0 * g_asymm);
  
alpha2 = -(1.0 - omega0) + 0.75 * (1.0 - omega0 * g_asymm);
  
 
lambda = sqrt (alpha1 * alpha1 - alpha2 * alpha2);
  
 
A =
    1.0 / (alpha2 / (alpha1 - lambda) * exp (lambda * dtau) -
	   alpha2 / (alpha1 + lambda) * exp (-lambda * dtau));
  
 
a11 = A * 2.0 * lambda / alpha2;
  
a12 = A * (exp (lambda * dtau) - exp (-lambda * dtau));
  
 
b = 0.5 - 0.75 * g_asymm * mu0;
  
alpha3 = -omega0 * b;
  
alpha4 = omega0 * (1 - b);
  
 
denom = (1.0 / mu0 / mu0 - lambda * lambda);
  
alpha5 = ((alpha1 - 1.0 / mu0) * alpha3 - alpha2 * alpha4) / denom;
  
alpha6 = (alpha2 * alpha3 - (alpha1 + 1.0 / mu0) * alpha4) / denom;
  
 
a33 = exp (-dtau / mu0);
  
 
a13 = alpha5 * (1.0 - (a11) * (a33)) - alpha6 * (a12);
  
a23 = -(a12) * alpha5 * (a33) + alpha6 * ((a33) - (a11));
  
 
*t = a11;
  
*r = a12;
  
*tdir = a33;
  
*rdir = a13 / mu0;
  
*sdir = a23 / mu0;

}


 
 
int
Define_Cloud (int layer, double tau_cloud_init, double *tau_cloud,
	      int cloud_layer)
{
  
for (int i = 0; i < layer; i++)
    {
      
tau_cloud[i] = 0;
    
} 
 
tau_cloud[cloud_layer] = tau_cloud_init;
  
return 0;

}


 
 
int
Magnus_formula (double *esat, double *T, int layer)
{
  
double Kconv = 273.15;
  
for (int il = 0; il < layer; il++)
    {
      
esat[il] =
	611.2 * exp ((17.43 * (T[il] - Kconv)) / (243.2 + (T[il] - Kconv)));
    
} 
return 0;

}


 
int
Calculat_init_relative_humidity (double *vmr, double *p, double *esat,
				 double *rel_hum, int layer)
{
  
 
for (int il = 0; il < layer; il++)
    {
      
rel_hum[il] = (vmr[il] * p[il]) / esat[il];
    
} 
return 0;

}


 
int
Rel_hum_To_vmr (double *vmr, double *esat, double *p, double *rel_hum,
		int layer)
{
  
for (int il = 0; il < layer; il++)
    {
      
vmr[il] = (rel_hum[il] * esat[il]) / p[il];
    
} 
return 0;

}


 
 
//--------------------------------------------- Functions End ------------------------------------------------------
  
 
 
// -------------------------------------------- Main function-------------------------------------------------------
  int
main ()
{
  
double *z = NULL;
  
double *p_lev = NULL;
  
double *T_lev = NULL;
  
double *h2ovmr_lev = NULL;
  
double *o3vmr_lev = NULL;
  
 
gnuplot_ctrl * g1;
  
g1 = gnuplot_init ();
  
 
    //parameter
  int level;
  
int k;
  
int ib;			//counting variable for bands
  int il;			// counting variable for layers
  int ny;			//size of ascii file
  int nwvlco2;			//dimension of wavelength vector Co2 ascii file
  int nwvlh2o;			//dimension of wavelength vector in H2O ascii file
  
double *wvl_lower_sw;
  
double *wvl_upper_sw;
  
double *qext_sw;
  
double *omega0_cloud_sw;
  
double *omega0_cloud_lw;
  
double *g_cloud_sw;
  
int nbands_sw;
  
 
double *wvl_lower_lw;
  
double *wvl_upper_lw;
  
double *qext_lw;
  
double *omega_lw;
  
double *g_cloud_lw;
  
int nbands_lw;
  
 
    /* read atmospheric profile of T and p */ 
  int status1 =
    read_5c_file ("fpda.atm", &z, &p_lev, &T_lev, &h2ovmr_lev, &o3vmr_lev,
		  &level);
  
int layer = level - 1;
  
 
 
    /* read in cloud parameters */ 
    // shortwave
  int status2 =
    read_5c_file ("./rrtm/cldprp/rrtm.sw.int", &wvl_lower_sw, &wvl_upper_sw,
		  &qext_sw, &omega0_cloud_sw, &g_cloud_sw, &nbands_sw);
  
    // longwave
  int status3 =
    read_5c_file ("./rrtm/cldprp/rrtm.lw.int", &wvl_lower_lw, &wvl_upper_lw,
		  &qext_lw, &omega0_cloud_lw, &g_cloud_lw, &nbands_lw);
  
 
 
    /* Convert z to meters */ 
    for (int i = 0; i < level; i++)
    {
      
z[i] = z[i] * 1000;
    
} 
 
    /* initialize volume mixing ratios */ 
  double co2vmr[layer], ch4vmr[layer], n2ovmr[layer], o2vmr[layer];
  
double cfc11vmr[layer], cfc12vmr[layer], cfc22vmr[layer], ccl4vmr[layer];
  
 
for (il = 0; il < layer; il++)
    {
      
co2vmr[il] = 400e-6;
      
ch4vmr[il] = 1.6e-6;
      
n2ovmr[il] = 320e-9;
      
o2vmr[il] = .209;
      
cfc11vmr[il] = 0;
      
cfc12vmr[il] = 0;
      
cfc22vmr[il] = 0;
      
ccl4vmr[il] = 0;
    
}
  
 
int nbands;		// number of bands
  double *band_lbound;		// [nbands]    
  double *band_ubound;		// [nbands]    
  double **dtau_mol_lw;		// [nbands][layer]       
  double **wgt_lw;		// [nbands][layer]
  double **dtau_mol_sw;		// [nbands][layer]       
  double **dtau_ray_sw;		// [nbands][layer]       
  double *wgt_sw;		// [nbands]  
  double tau_cloud[layer];
  
 
 
double T[layer];		//defining T[0], .. , T[9] for layer = 10
  double T_pot[layer];
  
double h2ovmr[layer];
  
double o3vmr[layer];
  
 
double delta_E[layer];
  
double delta_T[layer];
  
double p[layer];
  
double E_up[level];
  
double E_down[level];
  
double E_up_vis[level];
  
double E_down_vis[level];
  
double Edir[level];
  
double esat[layer];
  
double rel_hum[layer];
  
 
    //naturkonstanten
  double albedo = 0.3;
  
double E_solar = 340;
  
double T_ground = 288;
  
double mu0 = 0.25;		//easiest solution for sun elevation hight
  // cloud properties
  double tau_cloud_init = 1.4;	//tau of a cloud at layer #cloud_layer
  int cloud_layer = 19;		//cloud at layer #cloud_layer
  // ---------------------- Initial conditions ---------------------------------
  
    // Calculate T, h2ovmr, o3vmr, p from layer-vectors
    for (int i = 0; i < layer; i++)
    {
      
T[i] = 0.5 * (T_lev[i] + T_lev[i + 1]);
      
h2ovmr[i] = 0.5 * (h2ovmr_lev[i] + h2ovmr_lev[i + 1]) * 1e-6;
      
o3vmr[i] = 0.5 * (o3vmr_lev[i] + o3vmr_lev[i + 1]) * 1e-6;
      
p[i] = 0.5 * (p_lev[i] + p_lev[i + 1]) * 100;
    
} 
 
 
    /*Check Adiabatic initial profile
       printf("Check Adiabatic initial profile\n");
       for (int i = 0; i<layer; i++) {
       printf("%i\t%f\t%f\t%f\t%e\t%e\n", i, z[i], T[i],p[i],h2ovmr[i],o3vmr[i]);
       }
     */ 
    
    //Defines the cloud!
    Define_Cloud (layer, tau_cloud_init, tau_cloud, cloud_layer);
  
 
 
    // Calculate Initial relative humididty to keep constant during the timesteps by variation of T[i]  and rel_hum[init]--> vmr[i] goes into model
    Magnus_formula (esat, T, layer);
  
Calculat_init_relative_humidity (h2ovmr, p, esat, rel_hum, layer);
  
 
 
    // ---------------------- Start of Timeline ----------------------------------
    for (int TimeStep = 0; TimeStep < 100000; TimeStep++)
    {
      
 
 
	//calculates h2ovmr from changed T profilie and initial relative humidity file --> keep relative humidity constant
	Magnus_formula (esat, T, layer);
      
Rel_hum_To_vmr (h2ovmr, esat, p, rel_hum, layer);
      
 
	// ------------- Longwave Radiative Transfer ----------------------
	
	// Call C rrtm routine to get longwave absorption coefficients
	cfpda_rrtm_lw (layer, p_lev, T, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr,
		       o2vmr, 
cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr,
		       
&nbands, &band_lbound, &band_ubound, &wgt_lw,
		       &dtau_mol_lw);
      
 
Thermal_RT (T, layer, E_up, E_down, nbands, band_lbound, band_ubound,
		     wgt_lw, dtau_mol_lw, omega0_cloud_lw, tau_cloud);
      
 
 
 
	//Energy change in layer [i] is forwarded as temperature to layer [i]
	for (int i = 0; i < layer - 1; i++)
	{
	  
delta_E[i] = E_down[i] - E_up[i] + E_up[i + 1] - E_down[i + 1];
	  
delta_T[i] = delta_E[i] / 100;
	  
T[i] = T[i] + delta_T[i];
	
} 
 
 
delta_E[layer - 1] = E_down[layer - 1] - E_up[layer - 1];
      
delta_T[layer - 1] = delta_E[layer - 1] / 100;
      
T[layer - 1] = T[layer - 1] + delta_T[layer - 1];
      
 
 
	// ------------- Longwave Radiative Transfer END ----------------------
	
 
 
 
 
	// ------------- Shortwave Radiative Transfer ----------------------
	
 
	// Call C rrtm routine to get shortwave absorption coefficients
	cfpda_rrtm_sw (layer, p_lev, T, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr,
		       o2vmr, 
&nbands, &band_lbound, &band_ubound, &wgt_sw,
		       &dtau_mol_sw, &dtau_ray_sw);
      
 
 
 
 
Visual_RT (layer, E_up_vis, E_down_vis, Edir, nbands,
			  dtau_mol_sw, dtau_ray_sw, wgt_sw, mu0, albedo,
			  omega0_cloud_sw, tau_cloud, g_cloud_sw);
      
 
 
 
	//Energy change in layer [i] is forwarded as temperature to layer [i]
	for (int i = 0; i < layer - 1; i++)
	{
	  
delta_E[i] = E_down_vis[i] - E_up_vis[i] + E_up_vis[i + 1] - E_down_vis[i + 1] + Edir[i] - Edir[i + 1];	// Edir missing!!
	  delta_T[i] = delta_E[i] / 100;
	  
T[i] = T[i] + delta_T[i];
	
} 
 
 
 
 
delta_E[layer - 1] = E_down_vis[layer - 1] - E_up_vis[layer - 1] + Edir[layer - 1];	//!!??????? weg?? 
      delta_T[layer - 1] = delta_E[layer - 1] / 100;	// Edir missing
      T[layer - 1] = T[layer - 1] + delta_T[layer - 1];
      
 
	// ------------- Shortwave Radiative Transfer END----------------------
	
 
 
 
 
	// ----------------- Water Vapor Feedback Start-----------------
	
 
 
	// ----------------- Water Vapor Feedback End-----------------
	
 
 
 
 
 
 
 
 
	//-------------------- Advection ----------------------------
	
	//Calculate potential Temperature from Temperature profile
	To_T_pot (T, T_pot, layer, p);
      
 
	//Sort potential temperature profiles
	Sort_Layers (T_pot, layer);
      
 
	//Calculate back Temperature profiles form sorted potential Temperature
	T_pot_To_T (T, T_pot, layer, p);
      
 
	//Calculate the z-profile from Temperature and pressure
	HydrostaticEquation (layer, z, T, p);
      
 
	//-------------------- Advection ----------------------------
	
 
	//Ground Temperatur is set equal to last layer
	T_ground = T[layer - 1];
      
 
 
 
	//-------------------- Output ------------------------------
	
 
printf ("\nFlows @ t = %i\n", TimeStep);
      
for (int i = 0; i < layer; i++)
	{
	  
printf
	    ("i = %i\tz = %f\tp = %f\tT = %F\tE_down = %f\tE_up = %f \tE_down_vis = %f\tE_up_vis = %f\tEdir = %e\n",
	     i, z[i], p[i], T[i], E_down[i], E_up[i], E_down_vis[i],
	     E_up_vis[i], Edir[i]);
	
}
	
 
printf
	("i = %i\tz = 0\t\tp = %f\tT = %f\tE_down = %f\tE_up = %f\tE_down_vis = %f\tE_up_vis = %f\tEdir = %e\n",
	 layer, p_ground, T_ground, E_down[layer], E_up[layer],
	 E_down_vis[layer], E_up_vis[layer], Edir[layer]);
      
 
	//printf("\nTemperatur profile @ t = %i\n",TimeStep);
	for (int i = 0; i < layer; i++)
	{
	  
	    //print to Console
	    //printf("%u\t%u\t%f\t%f\t%f\t%f\t%f\n", TimeStep, i, z[i], delta_E[i], T[i], T_pot[i], T_ground);
	
} 
 
 
 
	//-------------------- Output End ------------------------------
	
 
 
 
	//-------------------- Live plotting ------------------------------
	
	/* plot every 10th time step */ 
	if (TimeStep % 10 == 0)
	{
	  
 
gnuplot_resetplot (g1);	/* start with new plot rather than plotting into exisiting one */
	  
gnuplot_setstyle (g1, "linespoints");	/* draw lines and points */
	  
gnuplot_set_xlabel (g1, "temperature [K]");	/* xaxis label */
	  
gnuplot_set_ylabel (g1, "altitude [km]");	/* yaxis label */
	  
 
 
	    /* plot temperature T as function of z and label with temperature */ 
	    gnuplot_plot_xy (g1, T, z, layer, "Temperature");
	
 
}
      
 
	//-------------------- Live plotting ENDE ------------------------------
    
 
 
}			// ---------------------- End of Timeline ----------------------------------
  
 
 
    /* close plot */ 
    gnuplot_close (g1);
  
 
printf ("program finished\n");
  
 
return 0;

 
}				// Main End