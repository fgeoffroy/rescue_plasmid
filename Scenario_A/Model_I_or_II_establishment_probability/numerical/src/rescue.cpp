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
  double sigma = i.sigma;
  double alpha = i.alpha;
  double N0 = i.N0;
  double s0 = i.s0;
  double tEst = i.tEst;

  // We integrate until the wildtype population is less than 1
  double tMax = log((s0 - N0 * alpha) / (N0 * (s0 - alpha))) / s0;



  //number of slices between 2*tEst and tMax
  int n = 5000;

  // Time array for Riemann sum
  double delta = (2 * tMax - tEst) / n;
  // We need to add time after tMax to compute the establishment probability of mutants which appears just before tMax
  // We add tMax, so that their establishment probability is computed on the same time duration as mutants appearing af t=0
  double tV[n];
  for (int i = 0; i < n; i++)
  {
    tV[i] = tEst + i * delta;  // Left Riemann sum
  }




  // Solving the system of ODE for each time point starting at t = tEst

  struct param_type_syst paramsSyst = {alpha, s0};
  gsl_odeiv2_system sys = {syst, NULL, 3, &paramsSyst};
  double epsabsSys = 1e-6;
  gsl_odeiv2_driver * d =
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk2,
                                  epsabsSys, epsabsSys, 0);

  double t = 0;
  double y[1] = { N0 };
  int status = gsl_odeiv2_driver_apply (d, &t, tEst, y);
  if (y[0] < 0) {
    y[0] = 0;
  }

  double NV[n];
  double muV[n];
  double rV[n];
  for (int i = 0; i < n; i++)
  {
    int status = gsl_odeiv2_driver_apply (d, &t, tV[i], y);
    y[0] < 0 ? NV[i] = 0 : NV[i] = y[0];
    muV[i] = 1 + c + alpha * NV[i];
    rV[i] = (1 + sM) * (1 - sigma) + beta * (1 - phi) * NV[i] - muV[i];
  }
  gsl_odeiv2_driver_free (d);




  // Establishment probabilities, mutation supply and rescue probability
  double sum1 = 0;
  for (int j = 0; j < n; j++)
  {
    double sum2 = 0;
    for (int k = 0; k <= j; k++)
    {
      sum2 += rV[k];
    }
    sum1 += muV[j] * exp(- delta * sum2);
  }

  double pEst0 = 1 / (1 + delta * sum1);

  std::ofstream log("resultsNumerical.txt", std::ios_base::app | std::ios_base::out);
  log << i.beta << " " << i.phi << " " << i.sM << " " << i.c << " " << i.sigma << " " << i.alpha << " " << i.N0 << " " << i.s0 << " " << i.tEst << " " << tMax << " " << n << " " << pEst0 << "\n";

  return 0;
}
