#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>
using namespace std;

// declaration and initialization of physical constants and model parameters
class Consts {
public: 
    static const double kappa; // adiabatic exponent [/]
    static const double c_air; // specific heat capacity [J/kg K]
    static const double g;     // gravity acceleration [m/s^2]
    static const double E_abs; // heating rate from surface [W/m^2]
    static const double sigma; // Stefan–Boltzmann constant [W/m^2 K^4]
    static const double eps;   // emissivity [/]
    
    static const int nlayer; // number of layers
    static const int nlevel; // number of levels
    static const double dt; // time step length [s]
    static const int n_steps; // number of timesteps [/]
    static const int output_steps; // intervall in which the model produces output [/]
    static const double cooling_rate; // prescribed cooling rate [K/s]
};

const double Consts::kappa = 2.0 / 7.0;
const double Consts::c_air = 1004;
const double Consts::g = 9.80665; 
const double Consts::E_abs = 235;
const double Consts::sigma = 5.670373e-8;
const double Consts::eps = 0.3;

const int Consts::nlayer = 25; 
const int Consts::nlevel = Consts::nlayer + 1; 
const double Consts::dt = 360.0; 
const int Consts::n_steps = 10000;
const int Consts::output_steps = 1000;
const double Consts::cooling_rate = - 1.475/ 86400.0; 

// first version of an output function, gets called for one timestep
void output(const float &time, const vector<double> &player,
            const vector<double> &T, const vector<double> &theta) {
  // print at what timestep the model is
  printf("Output at %2.1f hours \n", time);
  // print the pressure, temperature and potential temperature of each layer
  for (int i=0; i<player.size(); i++) {
    printf("%3d %6.1f %5.1f %5.1f\n",
           i, player[i], T[i], theta[i]);
  }
  return;
}

void t_to_theta(const vector<double> &T, vector<double> &theta, const vector<double> &conversion_factors) {
  transform(T.begin(), T.end(),
            conversion_factors.begin(), theta.begin(), multiplies<double>());
  return;
}

void theta_to_t(const vector<double> &theta, vector<double> &T, const vector<double> &conversion_factors) {
  transform(theta.begin(), theta.end(),
            conversion_factors.begin(), T.begin(), divides<double>());
  return;
}

// heating function - so far only surface heating
void thermodynamics(vector<double> &T, const double &dp) {
  // computation of surface heating
  T[T.size()-1]+= Consts::E_abs * Consts::dt * Consts::g / (Consts::c_air * dp * 100.0);
    
  // computation of thermal cooling at the ground (according to Stefan-Boltzmann law)
  T[T.size()-1]-= Consts::eps * pow(T[T.size()-1], 4) * Consts::sigma * Consts::dt * Consts::g / (Consts::c_air * dp * 100.0); 
    
  // prescribed cooling of atmosphere
  double delta_T = Consts::cooling_rate * Consts::dt;
  transform(T.begin(), T.end(), T.begin(), bind2nd(plus<double>(), delta_T));
  return;
}

int main() {
  double dp = 1000.0 / (double) Consts::nlayer;
  double dT = 100.0 / (double) Consts::nlayer;

  vector<double> player(Consts::nlayer); // vector of pressures for each layer
  vector<double> p(Consts::nlevel); // vector of pressures between the layers
  vector<double> T(Consts::nlayer); // vector of temperatures for each layer
  vector<double> theta(Consts::nlayer); // vector of pot. temperatures for each layer
  vector<double> conversion_factors(Consts::nlayer); // vector for the conversion factors between t and theta

  for (int i=0; i<Consts::nlevel; i++) {
    p[i] = dp * (double) i; // the pressure levels are spaced equally between 0 and 1000 hPa
  }

  for (int i=0; i<Consts::nlayer; i++) {
    T[i] = 180.0 + dT * (double) i; // just a first guess for the T-profile for each layer
    player[i] = (p[i]+p[i+1])/2.0; // computation of pressure between the levels
    conversion_factors[i] = pow(1000.0 / player[i], Consts::kappa); // computation of conversion factors
    theta[i] = T[i] * conversion_factors[i]; // computation of theta for each layer
  }

  /* end of initialization */

  // loop over time steps
  for (int i=0; i<=Consts::n_steps; i++) {
    // compute current time in hours from start
    float time = (float) i * Consts::dt / 360; // [hours]
    // calculate theta values from new T values
    t_to_theta(T, theta, conversion_factors);

    // sort theta to simulate a stabilizing mixing
    sort(theta.begin(), theta.end(), greater<double>());
    theta_to_t(theta, T, conversion_factors);

    // call output function every n=output_steps times
    if (i % Consts::output_steps == 0) {
      // call output function
      output(time, player, T, theta);
    }

    thermodynamics(T, dp);
  }

  return 0;
}
