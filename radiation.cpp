#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>
#include "cplkavg.h"
using namespace std;

// declaration and initialization of physical constants and model parameters
class Consts {
public: 
    static const double c_air; // specific heat capacity [J/kg K]
    static const double g;  // gravity acceleration [m/s^2]
    static const double T0; // sea level standard temperature [K]
    static const double M;  // molar mass of dry air [kg/mol]
    static const double R0; // universal gas constant [J/mol K]
    static const double h;  // Planck constant [J/s]
    static const double c;  // speed of light in vacuum [m/s]
    static const double kB; // Boltzmann constant [J/K]
    static const double sigma; // Stefan–Boltzmann constant [W/m^2 K^4]
    
    static const int nlayer; // number of layers
    static const int nlevel; // number of levels
    static const int nangle; // number of angles
    static const int nlamda; // number of wavelengths
    static const double tau_total; // total optical thickness of whole atmosphere
};

const double Consts::c_air = 1004.0;
const double Consts::g = 9.80665; 
const double Consts::T0 = 288.0;
const double Consts::M = 0.02896;
const double Consts::R0 = 8.3144;
const double Consts::h = 6.62607e-34;
const double Consts::c = 299792458;
const double Consts::kB = 1.380649e-23;
const double Consts::sigma = 5.670373e-8;

const int Consts::nlayer = 10; 
const int Consts::nlevel = Consts::nlayer + 1; 
const int Consts::nangle = 20; 
const int Consts::nlamda = 1000; 
const double Consts::tau_total = 1; 

// first version of an output function, gets called for one timestep
void output(const vector<double> &z, const vector<double> &p,
            const vector<double> &T, const vector<double> &E_down, const vector<double> &E_up) {
  // print the altitude, pressure and temperature of each level
  for (int i=0; i<Consts::nlevel; i++) {
    printf("%3d %6.3f %5.1f %5.2f %6.3f %6.3f\n",
           i, z[i], p[i], T[i], E_down[i], E_up[i]);
  }
  return;
}

// barometric formula (includes conversion to km)
double p_to_z(const double &p_level) {
  return Consts::c_air * Consts::T0 / Consts::g * (1 - pow(p_level / 1000.0, Consts::R0 / (Consts::c_air * Consts::M))) / 1000.0;
}

// Planck's law
double planck(const double &T_layer, const double &lamda_i) {
  return 2 * Consts::h * pow(Consts::c, 2) / (pow(lamda_i, 5) * (exp(Consts::h * Consts::c / (lamda_i * Consts::kB * T_layer)) - 1));
}

// Planck function for grey atmosphere integrated over the thermal spectral range
double grey_atmosphere(const double &T_layer) {
  return Consts::sigma * pow(T_layer, 4) / M_PI;
}

/*
// Planck function integrated over the thermal spectral range (sent by Pablo) 
double planck_integrated_infinity(const double &T_layer, const double &lamda){
  // compute powers of x, the dimensionless spectral coordinate
  double x = (Consts::h * Consts::c / Consts::kB) * lamda / T_layer;
  double x2 = pow(x, 2);
  double x3 = pow(x, 3);
  // decide how many terms of sum are needed
  double iterations = 2.0 + 20.0/x;
  iterations = (iterations<512) ? iterations : 512;
  int iter = int(iterations) ;
  // add up terms of sum
  double sum = 0;
  for (int n=1;  n<iter; n++){
    double dn = 1.0/n;
    sum += exp(-n * x) * (x3 + (3.0 * x2 + 6.0 * (x + dn) * dn) * dn) * dn;
  }
  return 2.0 * Consts::h * pow(Consts::c, 2) * pow(T_layer / (Consts::h * Consts::c / Consts::kB), 4) * sum;
}
double planck_integrated(const double &T_layer, const double &lamda0, const double &lamda1){
  double integral0 = planck_integrated_infinity(T_layer, lamda0);
  double integral1 = planck_integrated_infinity(T_layer, lamda1);
  return integral1 - integral0;
}
*/

// absorption coeficient 
double alpha (const double &tau, const double &mu){
  return 1.0 - exp(- tau / mu);
}

void thermal_radiative_transfer_monochromatic(vector<double> &T, vector<double> &E_down, vector<double> &E_up,
                                              vector<double> &tau, vector<double> &mu, const double &dmu, 
                                              const double &lamda0, const double &lamda1){
  for (int imu=0; imu<Consts::nangle; imu++) {
  // boundary conditions
  double L_down = 0.0; 
  double L_up = cplkavg(lamda0, lamda1, 288.0);
    for (int ilev=1; ilev<Consts::nlevel; ilev++) {
      L_down = (1 - alpha(tau[ilev - 1], mu[imu])) * L_down + alpha(tau[ilev - 1], mu[imu]) * cplkavg(lamda0, lamda1, T[ilev - 1]);
      E_down[ilev] += 2 * M_PI * L_down * mu[imu] * dmu;
    } 
    for (int ilev=Consts::nlevel-2; ilev >= 0; ilev--) {
      L_up = (1 - alpha(tau[ilev], mu[imu])) * L_up + alpha(tau[ilev], mu[imu]) * cplkavg(lamda0, lamda1, T[ilev]);
      E_up[ilev] += 2 * M_PI * L_up * mu[imu] * dmu;
    }
  } 
  return; 
}

int main() {

  double dp = 1000.0 / (double) Consts::nlayer;  
  double dmu = 1.0 / (double) Consts::nangle;   
  double dtau = Consts::tau_total / (double) Consts::nlayer;   
  double dlamda = 1000000.0 / (double) Consts::nlamda; // [nm]
  double lamda_interval[3][2] = {{1000,8000},{8000,12000},{12000,1000000}}; // [nm]

  vector<double> p(Consts::nlevel); // vector of pressures between the layers
  vector<double> z(Consts::nlevel); // vector of altitudes between the layers
  vector<double> Tlevel(Consts::nlevel); // vector of temperatures between the layers
  vector<double> T(Consts::nlayer); // vector of temperatures for each layer
  vector<double> tau(Consts::nlayer); // vector of optical thickness
  vector<double> mu(Consts::nangle); // vector of cosines of zenith angles, characterize direction of radiation
  vector<double> lamda(Consts::nlamda); // vector of wavelengths
  vector<double> E_down(Consts::nlayer); // vector of downgoing thermal irradiances for each layer
  vector<double> E_up(Consts::nlayer); // vector of upgoing thermal irradiances for each layer

    
  for (int i=0; i<Consts::nlevel; i++) {
    p[i] = dp * (double) i; // the pressure levels are spaced equally between 0 and 1000 hPa
    z[i] = p_to_z(p[i]); // altitude levels
    tau[i] = dtau;  // values of optical thickness for every layer are the same  
  }
  
  /*  
  for (int i=Consts::nlevel-1; i>=0; i--) {
    T[i] = Consts::T0 - 6.5 * z[i] ; // T-profile for each layer
  }
  */
    
  // given initial temperature profile
  Tlevel = {127.28, 187.09, 212.42, 229.22, 242.03, 252.48, 261.37, 269.13, 276.04, 282.29, 288.00};
    
  for (int i=0; i<Consts::nlayer; i++) {
    T[i] = (Tlevel[i] + Tlevel[i+1])/2.0; // T-profile for each layer
  }
    
  for (int i=0; i<Consts::nangle; i++) {
    mu[i] = dmu / 2.0 + dmu * (double) i; // angles are spaced equally between 0 and pi/2
                                          //mu[i] is the center of the i-interval
  }
    
  for (int i=0; i<Consts::nlamda; i++) {
    lamda[i] = 1000.0 + dlamda * (double) i; // [nm], wavelengths are spaced equally between 1 and 1000 microns 
  }

  for (int i=0; i<Consts::nlayer; i++) {
     E_down[i] = 0.0;
     E_up[i] = 0.0;
  } 
    
  /* end of initialization */
  
  for (int i=0; i<Consts::nlamda - 1; i++) {
    thermal_radiative_transfer_monochromatic(T, E_down, E_up, tau, mu, dmu, lamda[i], lamda[i+1]);
  } 
  
  output(z, p, T, E_down, E_up);
    
  return 0;   
}

