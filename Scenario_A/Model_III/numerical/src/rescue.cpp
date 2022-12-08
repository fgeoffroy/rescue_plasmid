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
  double c;
  double sigma;
  double beta;
};


int syst (double t, const double y[], double f[], void *params)
{
  (void)(t); /* avoid unused parameter warning */
  struct param_type_syst * p = (struct param_type_syst *)params;
  double alpha = p->alpha;
  double s0 = p->s0;
  double c = p->c;
  double sigma = p->sigma;
  double beta = p->beta;
  if (y[0] == 0) {
    f[0] = 0;
  } else {
    f[0] = y[0] * (s0 - alpha * (y[0] + y[1]) - beta * y[1]) + y[1] * (1 + s0) * sigma;
  }
  if (y[1] == 0) {
    f[1] = 0;
  } else {
    f[1] = y[1] * ((1 + s0) * (1 - sigma) - (1 + c) - alpha * (y[0] + y[1]) + beta * y[0]);
  }
  return GSL_SUCCESS;
}





int main(int argc, char const *argv[]) {

  rescueinput i=filetoinput(argv[1]);
  double beta = i.beta;
  double gamma = i.gamma;
  double sP = i.sP;
  double sM = i.sM;
  double c = i.c;
  double sigma = i.sigma;
  double alpha = i.alpha;
  double u = i.u;
  double s0 = i.s0;

  double betaM = gamma * beta;

  double NF0;
  double NP0;
  double beta_free = alpha * (c + sigma * (1 + sP)) / sP;
  if (beta <= beta_free) {
    NF0 = sP / alpha;
    NP0 = 0;
  } else {
    if (sigma == 0) {
      double beta_wild = c * alpha / (sP - c);
      if (beta >= beta_wild) {
        NF0 = 0;
        NP0 = (sP - c) / alpha;
      } else {
        NF0 = (alpha * c - beta * (sP - c)) / (beta * beta);
        NP0 = (beta * sP - alpha * c) / (beta * beta);
      }
    } else {
      NF0 = (alpha*c + sqrt(c*(pow(alpha + beta, 2)*c + 4*alpha*beta*sigma) - 2*beta*c*(alpha + beta - 2*alpha*sigma)*sP + pow(beta*sP,2)) + beta*(c - sP + 2*sigma*(1 + sP)))/(2.*pow(beta,2));
      NP0 = -(pow(alpha,2)*c + alpha*(-(beta*sP) + 2*beta*sigma*(1 + sP) + sqrt(pow(alpha*c,2) + pow(beta*(c-sP),2) + 2*alpha*beta*c*(c - sP + 2*sigma*(1 + sP)))) - beta*(beta*(-c + sP) + sqrt(pow(alpha*c,2) + pow(beta*(c-sP),2) + 2*alpha*beta*c*(c - sP + 2*sigma*(1 + sP)))))/(2.*alpha*pow(beta,2));
    }
  }

  double N0 = NF0 + NP0;

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
  double NFV[ntot];
  double NPV[ntot];
  double muV[ntot];
  double rV[ntot];
  struct param_type_syst paramsSyst = {alpha, s0, c, sigma, beta};
  gsl_odeiv2_system sys = {syst, NULL, 2, &paramsSyst};
  double epsabsSys = 1e-6;
  gsl_odeiv2_driver * d =
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk2,
                                  epsabsSys, epsabsSys, 0);

  double t = 0;
  double y[2] = { NF0, NP0 };
  for (int i = 0; i < ntot; i++)
  {
    int status = gsl_odeiv2_driver_apply (d, &t, tV[i], y);
    y[0] < 0 ? NFV[i] = 0 : NFV[i] = y[0];
    y[1] < 0 ? NPV[i] = 0 : NPV[i] = y[1];
    muV[i] = 1 + c + alpha * (NFV[i] + NPV[i]);
    rV[i] = (1 + sM) * (1 - sigma) + betaM * NFV[i] - muV[i];
  }
  gsl_odeiv2_driver_free (d);




  // Establishment probabilities, mutation supply and rescue probability
  double pEstV[n];
  double mutV[n];
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
    sum1 += u * (1 + s0) * NPV[i] * ( 1 / (1 + delta * sum2) );
  }
  double pRes = 1 - exp(- delta * sum1);

  std::ofstream log("resultsNumerical.txt", std::ios_base::app | std::ios_base::out);
  log << i.beta << " " << i.gamma << " " << i.sP << " " << i.sM << " " << i.c << " " << i.sigma << " " << i.alpha << " " << i.u << " " << i.s0 << " " << tMax << " " << n << " " << pRes << "\n";

  return 0;
}
