#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>
using namespace std;

void output(const float &time, const vector<double> &player,
            const vector<double> &T, const vector<double> &theta) {
  printf("Output at %2.1f hours \n", time);
  for (int i=0; i<player.size(); i++) {
    printf("%3d %6.1f %5.1f %5.1f\n",
           i, player[i], T[i], theta[i]);
  }
  return;
}

void heating(vector<double> &theta, const double &dp, const double &heatrate,
             const double &dt, const double &g, const double &c_air,
             const vector<double> &t_to_theta) {
  double delta_T = heatrate * dt * g / (c_air * dp);
  theta[theta.size()-1] += delta_T * t_to_theta[theta.size()-1];
  return;
}

int main() {
  int nlayer = 10;
  double kappa = 2.0 / 7.0; // adiabatic exponent [/]
  double c_air = 1004; // specific heat capacity [J/kg K]
  double g = 9.81; // gravity acceleration [m/s^2]
  double dt = 360.0; // time step length [s]
  double n_steps = 4; // experiment lenth [s]
  int output_steps = 1;
  double heatrate = 2; // heating rate from surface [W/m^2]
  
  double dp = 1000.0 / (double) nlayer;
  double dT = 100.0 / (double) nlayer;
  
  int nlevel = nlayer + 1;
  vector<double> player(nlayer);
  vector<double> p(nlevel);
  vector<double> T(nlayer);
  vector<double> theta(nlayer);
  vector<double> t_to_theta(nlayer);

  for (int i=0; i<nlevel; i++) {
    p[i] = dp * (double) i;
  }
  
  for (int i=0; i<nlayer; i++) {
    T[i] = 180.0 + dT * (double) i;
    player[i] = (p[i]+p[i+1])/2.0;
    t_to_theta[i] = pow(1000.0 / player[i], kappa);
    theta[i] = T[i] * t_to_theta[i];
    // printf("%3d %6.1f %5.1f %5.1f\n", i, player[i], T[i], theta[i]);
  }
  

  /* end of initialization */

// loop over time steps
for (int i=0; i<=n_steps; i++) {
  float time = (float) i * dt / 360; // [hours]

  // Atmosphere stabilizes
  sort(theta.begin(), theta.end(), greater<double>());

  if (i % output_steps == 0) {
    // calculate temperature
    transform(theta.begin(), theta.end(), t_to_theta.begin(), T.begin(), divides<double>());
    output(time, player, T, theta);
  }

  heating(theta, dp, heatrate, dt, g, c_air, t_to_theta);
}
  
return 0;
}