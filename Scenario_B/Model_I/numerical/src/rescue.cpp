#include <array>
#include <gsl/gsl_rng.h>
#include <algorithm>
#include <iterator>
#include <math.h>
#include <gsl/gsl_odeiv2.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <cstdio>
#include "include/fileio.h"
using namespace std;







struct param_type_syst {
  double alpha;
  double s0;
};


int syst (double t, const double y[], double f[], void *params)
{
  (void)(t); /* avoid unused parameter warning */
  struct param_type_syst * p = (struct param_type_syst *)params;
  double alpha  = p->alpha;
  double s0  = p->s0;
  if (y[0] == 0) {
    f[0] = 0;
  } else {
    f[0] = y[0] * (s0 - alpha * y[0]);
  }
  return GSL_SUCCESS;
}






int main(int argc, char const *argv[]) {

  rescueinput i=filetoinput(argv[1]);
  double beta = i.beta;
  double phi = i.phi;
  double sM = i.sM;
  double c = i.c;
  double alpha = i.alpha;
  double u = i.u;
  double N0 = i.N0;
  double s0 = i.s0;

  // We integrate until the wildtype population is less than 1
  double tMax = log((s0 - N0 * alpha) / (N0 * (s0 - alpha))) / s0;
  //number of slices between 0 and tMax
  int n = 5000;

  // Time array for Riemann sum
  double delta = tMax / n;
  // We need to add time after tMax to compute the establishment probability of mutants which appears just before tMax
  // We add tMax, so that their establishment probability is computed on the same time duration as mutants appearing af t=0
  int n2 = n;
  int ntot = n + n2;
  double tV[ntot];
  for (int i = 0; i < ntot; i++)
  {
    tV[i] = i * delta;  // Left Riemann sum
  }


  // Solving the system of ODE for each time point
  double NV[ntot];
  double muV[ntot];
  double rV[ntot];
  struct param_type_syst paramsSyst = {alpha, s0};
  gsl_odeiv2_system sys = {syst, NULL, 3, &paramsSyst};
  double epsabsSys = 1e-6;
  gsl_odeiv2_driver * d =
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk2,
                                  epsabsSys, epsabsSys, 0);

  double t = 0;
  double y[1] = { N0 };
  for (int i = 0; i < ntot; i++)
  {
    int status = gsl_odeiv2_driver_apply (d, &t, tV[i], y);
    y[0] < 0 ? NV[i] = 0 : NV[i] = y[0];
    muV[i] = 1 + c + alpha * NV[i];
    rV[i] = 1 + sM + beta * (1 - phi) * NV[i] - muV[i];
  }
  gsl_odeiv2_driver_free (d);




  // Establishment probabilities, mutation supply and rescue probability
  double sum1 = 0;
  for (int i = 0; i < n; i++)
  {
    double sum2 = 0;
    for (int j = i; j < ntot; j++)
    {
      double sum3 = 0;
      for (int k = i; k <= j; k++)
      {
        sum3 += rV[k];
      }
      sum2 += muV[j] * exp(- delta * sum3);
    }
    sum1 += u * (1 - phi) * NV[i] * ( 1 / (1 + delta * sum2) );
  }
  double pRes = 1 - exp(- delta * sum1);

  std::ofstream log("resultsNumerical.txt", std::ios_base::app | std::ios_base::out);
  log << i.beta << " " << i.phi << " " << i.sM << " " << i.c << " " << i.alpha << " " << i.u << " " << i.N0 << " " << i.s0 << " " << tMax << " " << n << " " << pRes << "\n";

  return 0;
}
