/* cplkavg was extracted from cdisort.c which is part of libRadtran; */
/* c_planck_func () is identical to c_planck_func1() in cdisort.c;   */
/* Bernhard Mayer, 23.11.2020                                        */
/* Addapted for C++, Tatsiana Bardachova, 30.11.2020                 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <string>
#include <ctype.h>
#include <limits.h>
#include <float.h>


#define DS_WARNING 0
#define DS_ERROR 1


/*============================= c_errmsg() ===============================*/

/*
 * Print out a warning or error message;  abort if type == DS_ERROR
 */

#define MAX_WARNINGS 100

void c_errmsg (std::string messag, int type) {
  static int warning_limit = 0, num_warnings = 0;

  if (type == DS_ERROR) {
    fprintf (stderr, "\n ******* ERROR >>>>>>  %s\n", messag.c_str());
    exit (1);
  }

  if (warning_limit)
    return;

  if (++num_warnings <= MAX_WARNINGS) {
    fprintf (stderr, "\n ******* WARNING >>>>>>  %s\n", messag.c_str());
  } else {
    fprintf (stderr, "\n\n >>>>>>  TOO MANY WARNING MESSAGES --  ','They will no longer be printed  <<<<<<<\n\n");
    warning_limit = 1;
  }

  return;
}

#undef MAX_WARNINGS

/*============================= end of c_errmsg() ========================*/





/*============================= c_planck_func() ========================*/

/*
  Computes Planck function integrated between two wavelengths
  INPUT :  wvllo  : Lower wavelength (nm) of spectral interval
           wvlhi  : Upper wavelength (nm)
           t      : Temperature (K)
  OUTPUT : ans    : Integrated Planck function ( Watts/sq m ) = Integral (WNUMLO to WNUMHI) of
                    2h c*c nu*nu*nu / ( exp(hc nu/kT) - 1) (where h=Plancks constant, c=speed of
                    light, nu=wavenumber, T=temperature, and k = Boltzmann constant)
  Reference : Specifications of the Physical World: New Value of the Fundamental Constants,
                Dimensions/N.B.S., Jan. 1974
  Method :  For wnumlo close to wnumhi, a Simpson-rule quadrature is done to avoid
            ill-conditioning; otherwise
            (1)  For wnumlo or wnumhi small, integral(0 to WNUMLO/HI) is calculated by expanding
                 the integrand in a power series and integrating term by term;
            (2)  Otherwise, integral (wnumlo/hi to infinity) is calculated by expanding the
                 denominator of the integrand in powers of the exponential and integrating
                 term by term.
  Accuracy :  At least 6 significant digits, assuming the physical constants are infinitely accurate
  ERRORS WHICH ARE NOT TRAPPED:
      * Power or exponential series may underflow, giving no significant digits.  
        This may or may not be of concern, depending on the application.
      * Simpson-rule special case is skipped when denominator of integrand will cause overflow.
        In that case the normal procedure is used, which may be inaccurate if the wavenumber limits
        (wnumlo, wnumhi) are close together.
  LOCAL VARIABLES
        a1,2,... :  Power series coefficients
        c2       :  h * c / k, in units cm*K (h = Plancks constant,
                      c = speed of light, k = Boltzmann constant)
        D(I)     :  Exponential series expansion of integral of
                       Planck function from wnumlo (i=1) or wnumhi
                       (i=2) to infinity
        ex       :  exp( - V(I) )
        exm      :  ex**m
        mmax     :  No. of terms to take in exponential series
        mv       :  Multiples of V(I)
        P(I)     :  Power series expansion of integral of
                       Planck function from zero to WNUMLO (I=1) or
                       WNUMHI (I=2)
        sigma    :  Stefan-Boltzmann constant (W/m**2/K**4)
        sigdpi   :  sigma/pi
        smallv   :  Number of times the power series is used (0,1,2)
        V(i)     :  c2 * (wnumlo(i=1) or wnumhi(i=2))/temperature
        vcut     :  Power-series cutoff point
        vcp      :  Exponential series cutoff points
        vmax     :  Largest allowable argument of EXP function
   Called by- c_disort
   Calls- c_errmsg
 ----------------------------------------------------------------------*/

#define PLKF(x)                                                                                                                    \
  ({                                                                                                                               \
    const double _x = (x);                                                                                                         \
    _x* _x* _x / (exp (_x) - 1.);                                                                                                  \
  })
#define A1 (1. / 3.)
#define A2 (-1. / 8.)
#define A3 (1. / 60.)
#define A4 (-1. / 5040.)
#define A5 (1. / 272160.)
#define A6 (-1. / 13305600.)
#define C2 (1.438786)
#define SIGMA (5.67032E-8)
#define VCUT (1.5)

double cplkavg (double wvllo, double wvlhi, double t) {
  register int  i, k, m, mmax, n, smallv;
  int           converged;
  static int    initialized = 0;
  const double  vcp[7]      = {10.25, 5.7, 3.9, 2.9, 2.3, 1.9, 0.0};
  double        del, ex, exm, hh, mv, oldval, val, val0, vsq, d[2], p[2], v[2], ans;
  static double vmax, sigdpi, conc;

  if (!initialized) {
    vmax   = log (DBL_MAX);
    sigdpi = SIGMA / M_PI;
    conc   = 15. / pow (M_PI, 4.);

    initialized = 1;
  }

  /* convert from wavelength to wavenumber */
  double wnumhi = 1.0E7 / wvllo;
  double wnumlo = 1.0E7 / wvlhi;
  
  if (t < 0. || wnumhi <= wnumlo || wnumlo < 0.) {
    c_errmsg ("planck_func1--temperature or wavenums. wrong", DS_ERROR);
  }

  if (t < 1.e-4) {
    return 0.;
  }

  v[0] = C2 * wnumlo / t;
  v[1] = C2 * wnumhi / t;

  if (v[0] > DBL_EPSILON && v[1] < vmax && (wnumhi - wnumlo) / wnumhi < 1.e-2) {
    /*
     * Wavenumbers are very close.  Get integral by iterating Simpson rule to convergence.
     */
    hh        = v[1] - v[0];
    oldval    = 0.;
    val0      = PLKF (v[0]) + PLKF (v[1]);
    converged = 0;
    for (n = 1; n <= 10; n++) {
      del = hh / (2 * n);
      val = val0;
      for (k = 1; k <= 2 * n - 1; k++) {
        val += (double)(2 * (1 + k % 2)) * PLKF (v[0] + (double)k * del);
      }
      val *= del * A1;
      if (fabs ((val - oldval) / val) <= 1.e-6) {
        /* convergence */
        converged = 1;
        break;
      }
      oldval = val;
    }
    if (!converged) {
      c_errmsg ("planck_func1--Simpson rule didn't converge", DS_WARNING);
    }

    return sigdpi * pow (t, 4.0) * conc * val;
  }

  /*
   * General case
   */
  smallv = 0;
  for (i = 1; i <= 2; i++) {
    if (v[i - 1] < VCUT) {
      /*
       * Use power series
       */
      smallv++;
      vsq      = (v[i - 1]) * (v[i - 1]);
      p[i - 1] = conc * vsq * v[i - 1] * (A1 + v[i - 1] * (A2 + v[i - 1] * (A3 + vsq * (A4 + vsq * (A5 + vsq * A6)))));
    } else {
      /*
       * Use exponential series
       *
       * Find upper limit of series
       */
      mmax = 1;
      while (v[i - 1] < vcp[mmax - 1]) {
        mmax++;
      }

      ex       = exp (-v[i - 1]);
      exm      = 1.;
      d[i - 1] = 0.;
      for (m = 1; m <= mmax; m++) {
        mv  = (double)m * v[i - 1];
        exm = ex * exm;
        d[i - 1] += exm * (6. + mv * (6. + mv * (3. + mv))) / (m * m * m * m);
      }
      d[i - 1] *= conc;
    }
  }
  /*
   * Handle ill-conditioning
   */
  if (smallv == 2) {
    /*
     * wnumlo and wnumhi both small
     */
    ans = p[1] - p[0];
  } else if (smallv == 1) {
    /*
     * wnumlo small, wnumhi large
     */
    ans = 1. - p[0] - d[1];
  } else {
    /*
     * wnumlo and wnumhi both large
     */
    ans = d[0] - d[1];
  }
  ans *= sigdpi * pow (t, 4.0);
  if (ans == 0.) {
    c_errmsg ("planck_func1--returns zero; possible underflow", DS_WARNING);
  }

  return ans;
}

#undef A1
#undef A2
#undef A3
#undef A4
#undef A5
#undef A6
#undef C2
#undef SIGMA
#undef VCUT
#undef PLKF

/*============================= end of c_planck_func() =================*/
