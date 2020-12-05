#include <stdio.h>
#include <math.h>

#include "cplkavg.h"

int main () {
  double T=288.0;
  double inc = 1.01;

  // print integrated Planck function times pi for a number of intervals
  // between 0 1000 nm and 1000 micron; sum should equal sigma * T**4

  for (double wvl=1000;wvl<1000000; wvl*=inc) {
    fprintf (stdout, "%10.2f %10.2f  %e\n", wvl, wvl*inc,
	     M_PI*cplkavg (wvl, wvl*inc, T));
  }

  return 0;
}

