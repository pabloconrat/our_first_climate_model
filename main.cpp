#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>
using namespace std;

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
void thermodynamics(vector<double> &T, const double &dp, const double &E_abs,
             const double &dt, const double &g, const double &c_air, const double cooling_rate,
             const vector<double> &conversion_factors) {
  // computation of surface heating
  T[T.size()-1]+= E_abs * dt * g / (c_air * dp * 100.0);
  
  // prescribed cooling of atmosphere
  double delta_T = cooling_rate * dt;
  transform(T.begin(), T.end(), T.begin(), bind2nd(plus<double>(), delta_T));
  return;
}

// declaration and initialization of physical constants 
class consts {
public: 
    static const double kappa; // adiabatic exponent [/]
    static const double c_air; // specific heat capacity [J/kg K]
    static const double g;     // gravity acceleration [m/s^2]
    static const double E_abs; // heating rate from surface [W/m^2]
};

const double consts::kappa = 2.0 / 7.0;
const double consts::c_air = 1004;
const double consts::g = 9.80665; 
const double consts::E_abs = 235;


int main() {
    
  int nlayer = 25; // number of layers
  double dt = 360.0; // time step length [s]
  double n_steps = 1000; // number of timesteps [/]
  int output_steps = 100; // intervall in which the model produces output [/]
  double cooling_rate = - 3.0 / 86400.0; // prescribed cooling rate [K/s]

  double dp = 1000.0 / (double) nlayer;
  double dT = 100.0 / (double) nlayer;

  int nlevel = nlayer + 1;
  vector<double> player(nlayer); // vector of pressures for each layer
  vector<double> p(nlevel); // vector of pressures between the layers
  vector<double> T(nlayer); // vector of temperatures for each layer
  vector<double> theta(nlayer); // vector of pot. temperatures for each layer
  vector<double> conversion_factors(nlayer); // vector for the conversion factors between t and theta

  for (int i=0; i<nlevel; i++) {
    p[i] = dp * (double) i; // the pressure levels are spaced equally between 0 and 1000 hPa
  }

  for (int i=0; i<nlayer; i++) {
    T[i] = 180.0 + dT * (double) i; // just a first guess for the T-profile for each layer
    player[i] = (p[i]+p[i+1])/2.0; // computation of pressure between the levels
    conversion_factors[i] = pow(1000.0 / player[i], consts::kappa); // computation of conversion factors
    theta[i] = T[i] * conversion_factors[i]; // computation of theta for each layer
  }

  /* end of initialization */

  // loop over time steps
  for (int i=0; i<=n_steps; i++) {
    // compute current time in hours from start
    float time = (float) i * dt / 360; // [hours]
    // calculate theta values from new T values
    t_to_theta(T, theta, conversion_factors);

    // sort theta to simulate a stabilizing mixing
    sort(theta.begin(), theta.end(), greater<double>());

    // call output function every n=output_steps times
    if (i % output_steps == 0) {
      // call output function
      output(time, player, T, theta);
    }

    thermodynamics(T, dp, consts::E_abs, dt, consts::g, consts::c_air, cooling_rate, conversion_factors);
  }

  return 0;
}
