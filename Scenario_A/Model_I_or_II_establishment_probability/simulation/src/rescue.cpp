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



// main function for rescue simulation
bool rescue(double beta, double phi, double sM, double c, double sigma, double alpha,
  unsigned int N0, double s0, unsigned int seed) {

  // DECLARE+INITIALIZE VARIABLES
  // initialize variables for gillespie algorithm (see Barnes(2015) eq. 7.2)
  double t=0.0; // set time to zero
  double r1,r2,r3;   // random numbers
  // auxiliary variables
  int r=0;  // next reaction index
  int i; // Individual Class with i mutated plasmids
  double sum_a_i; // aux. variable for gillespie algorithm
  double a0; // sum of a_i values
  // population size at time t=0
  unsigned long long Nvec[2];
  Nvec[0]=N0;
  Nvec[1]=1;  //We start with one mutant to study its establishment probability
  int N=Nvec[0]+Nvec[1]; // pop. size
  // propensity of birth, death and transfer
  double a[5]; // declare array for propensity values.
  // declare random distributions
  const gsl_rng_type * T;
  gsl_rng * rng1;
  T = gsl_rng_default;
  rng1 = gsl_rng_alloc (T);
  gsl_rng_set(rng1,seed);

  // SIMULATION

  // For the case that rescue plasmid cells have a negative vertical fitness, return 0.
  if ( (1 + c) >= (1 + sM) * (1 - sigma) ) {
    return 0;
  }

  while ((Nvec[1]!=0) && Nvec[1]<= log (.01) / log((1 + c) / ((1 + sM) * (1 - sigma))) ) {
    // refresh propensity values a[]
    /** birth **/
    a[0] = Nvec[0] * (1 + s0);
    a[1] = Nvec[1] * (1 + sM) * (1 - sigma);
    /** death **/
    a[2] = Nvec[0] * (1 + alpha * N);
    a[3] = Nvec[1] * (1 + c + alpha * N);
    /** horizontal transfer **/
    a[4] = beta * (1 - phi) * Nvec[0] * Nvec[1];
    a0=0.0; for (auto& ai : a) a0 += ai; // sum of a[] of all processes
    // generate random numbers from uniform distribution
    r1=gsl_rng_uniform_pos(rng1);
    r2=gsl_rng_uniform_pos(rng1);
    // identify next reaction
    sum_a_i=0.0;r=0;
    do {sum_a_i += a[r];
      if (r2*a0<=sum_a_i) break;
      r++;}
      while(r<5);
    // increase current time
    t+=1.0/a0*log(1.0/r1);
    // evaluate reaction r
    if (r<2) { // birth
      if (r == 0) {
        Nvec[0]++;
      } else {
        Nvec[1]++;
      }
      N++;
    }
    else if (r < 4) { // death
      if (r == 2) {
        Nvec[0]--;
      } else {
        Nvec[1]--;
      }
      N--;
    } else { // transfer of a mutant plasmid
      Nvec[0]--;
      Nvec[1]++;
    }

  }
  // return : 0 <-> population extinction, 1 <-> population rescue
  return (Nvec[1]!=0);
}




int main(int argc, char const *argv[]) {

  rescueinput i=filetoinput(argv[1]);
  int NSim = 10000;


  struct param_type_syst paramsSyst = {i.alpha, i.s0};
  gsl_odeiv2_system sys = {syst, NULL, 3, &paramsSyst};
  double epsabsSys = 1e-6;
  gsl_odeiv2_driver * d =
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk2,
                                  epsabsSys, epsabsSys, 0);

  double t = 0;
  double y[1] = { i.N0 };
  int status = gsl_odeiv2_driver_apply (d, &t, i.tEst, y);
  if (y[0] < 0) {
    y[0] = 0;
  }
  gsl_odeiv2_driver_free (d);


  double result=0.;
  for (size_t k = 0; k < NSim; k++) {
    bool b=rescue(i.beta, i.phi, i.sM, i.c, i.sigma, i.alpha, (int)round(y[0]), i.s0, k);
    result+= b;
  }
  result/=NSim;

  std::ofstream log("resultsSimulation.txt", std::ios_base::app | std::ios_base::out);
  log << i.beta << " " << i.phi << " " << i.sM << " " << i.c << " " << i.sigma << " " << i.alpha << " " << i.N0 << " " << i.s0 << " " << i.tEst << " " << NSim << " " << result << "\n";

  return 0;
}
