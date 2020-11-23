#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>
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
const int Consts::nangle = 10; 
const int Consts::nlamda = 1000; 

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

// Planck function for 'window' atmosphere integrated over the thermal spectral range
double window_atmosphere(const double &T_layer, const double &lamda0, const double &lamda1, const double &dlamda) { 
  double integral = 0.0;
  int nlamda = (int) ceil (lamda1 - lamda0) / dlamda;
  for (int i=0; i<nlamda; i++) {
    double lamda = lamda0 + dlamda * (double) i;
    integral += planck(T_layer, lamda) * dlamda;
  }
  return integral;
}

// absorption coeficient 
double alpha (const double &tau, const double &mu){
  return 1.0 - exp(-tau / mu);
}

void radiative_transfer(vector<double> &B, vector<double> &E_down, vector<double> &E_up, 
                        vector<double> &tau, vector<double> &mu, const double &dmu){
  for (int imu=0; imu<Consts::nangle; imu++) {
  // boundary conditions
  double L_down = 0.0; 
  double L_up = B[B.size()-1]; // ???
    for (int ilev=1; ilev<Consts::nlevel; ilev++) {
      L_down = (1 - alpha(tau[ilev - 1], mu[imu])) * L_down + alpha(tau[ilev - 1], mu[imu]) * B[ilev - 1];
      E_down[ilev] += 2 * M_PI * L_down * mu[imu] * dmu;
    } 
    for (int ilev=Consts::nlevel-2; ilev >= 0; ilev--) {
      L_up = (1 - alpha(tau[ilev], mu[imu]))*L_up + alpha(tau[ilev], mu[imu])*B[ilev];
      E_up[ilev] += 2 * M_PI * L_up * mu[imu] * dmu;;
    }
  } 
  return; 
}

int main() {

  double dp = 1000.0 / (double) Consts::nlayer;  
  double dmu = 1.0 / (double) Consts::nangle;    
  double dlamda = 1000e-6 / (double) Consts::nlamda; 
  double lamda_interval[3][2] = {{7e-9,8e-6},{8e-6,12e-6},{12e-6,1000e-6}};

  vector<double> p(Consts::nlevel); // vector of pressures between the layers
  vector<double> z(Consts::nlevel); // vector of altitudes between the layers
  vector<double> T(Consts::nlayer); // vector of temperatures for each layer
  vector<double> tau(Consts::nlayer); // vector of optical thickness
  vector<double> mu(Consts::nangle); // vector of cosines of zenith angles, characterize direction of radiation
  vector<double> lamda(Consts::nlamda); // vector of wavelengths
  vector<double> B(Consts::nlayer); // vector of spectral radiances for each layer
  vector<double> E_down(Consts::nlayer); // vector of downgoing thermal irradiances for each layer
  vector<double> E_up(Consts::nlayer); // vector of upgoing thermal irradiances for each layer

    
  for (int i=0; i<Consts::nlevel; i++) {
    p[i] = dp * (double) i; // the pressure levels are spaced equally between 0 and 1000 hPa
    z[i] = p_to_z(p[i]); // altitude levels
    tau[i] = (double) i + 1.0;  // ???
  }
    
  for (int i=Consts::nlevel-1; i>=0; i--) {
    T[i] = Consts::T0 - 6.5 * z[i] ; // T-profile for each layer
  }
    
  for (int i=0; i<Consts::nangle; i++) {
    mu[i] = dmu / 2.0 + dmu * (double) i; // angles are spaced equally between 0 and pi/2
                                          //mu[i] is the center of the i-interval
  }
    
  for (int i=0; i<Consts::nlamda; i++) {
    lamda[i] = 1.0 + dlamda * (double) i; // wavelengths are spaced equally between 1 and 15 microns
  }

  for (int i=0; i<Consts::nlayer; i++) {
     E_down[i] = 0.0;
     E_up[i] = 0.0;
  } 
    
  /* end of initialization */
  
  // grey atmosphere
  for (int i=0; i<Consts::nlayer; i++) {
    B[i] = grey_atmosphere(T[i]);
  } 
   
  radiative_transfer(B, E_down, E_up, tau, mu, dmu);
    
  output(z, p, T, E_down, E_up);
    
  return 0;   
}

