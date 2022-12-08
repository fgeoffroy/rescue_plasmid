#include <array>
#include <gsl/gsl_rng.h>
#include <algorithm>
#include <iterator>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <cstdio>
#include "include/fileio.h"
using namespace std;

// main function for rescue simulation
bool rescue(double beta, double phi0, double gamma, double sM, double c,
  double sigma, double alpha, double u, unsigned int N0, double s0, unsigned int seed) {

  // DECLARE+INITIALIZE VARIABLES
  double betaM = gamma * beta;
  // initialize variables for gillespie algorithm (see Barnes(2015) eq. 7.2)
  double t=0.0; // set time to zero
  double r1,r2,r3;   // random numbers
  // auxiliary variables
  int r=0;  // next reaction index
  int i; // Individual Class with i mutated plasmids
  double sum_a_i; // aux. variable for gillespie algorithm
  double a0; // sum of a_i values
  // population size at the t=0
  unsigned long long Nvec[3];
  Nvec[0]=round((1 - phi0) * N0);
  Nvec[1]=N0-Nvec[0];
  Nvec[2]=0;
  int N=Nvec[0]+Nvec[1]+Nvec[2]; // pop. size
  // propensity of birth, death and transfer
  double a[9]; // declare array for propensity values.
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

  while ((N!=0) && Nvec[2]<= log (.01) / log((1 + c) / ((1 + sM) * (1 - sigma))) ) {
    // refresh propensity values a[]
    /** birth **/
    a[0] = Nvec[0] * (1 + s0);
    a[1] = Nvec[1] * (1 + s0);
    a[2] = Nvec[2] * (1 + sM);
    /** death **/
    a[3] = Nvec[0] * (1 + alpha * N);
    a[4] = Nvec[1] * (1 + c + alpha * N);
    a[5] = Nvec[2] * (1 + c + alpha * N);
    /** horizontal transfer **/
    a[6] = beta * Nvec[0] * Nvec[1];
    a[7] = betaM * Nvec[0] * Nvec[2];
    /** Immigration event of an incompatible resistant plasmid via HT (we assume no effect of beta here) **/
    a[8] = u * Nvec[0];
    a0=0.0; for (auto& ai : a) a0 += ai; // sum of a[] of all processes
    // generate random numbers from uniform distribution
    r1=gsl_rng_uniform_pos(rng1);
    r2=gsl_rng_uniform_pos(rng1);
    // identify next reaction
    sum_a_i=0.0;r=0;
    do {sum_a_i += a[r];
      if (r2*a0<=sum_a_i) break;
      r++;}
      while(r<9);
    // increase current time
    t+=1.0/a0*log(1.0/r1);
    // evaluate reaction r
    if (r<3) { // birth
      if (r == 0) {
        Nvec[0]++;
      } else {
        r3=gsl_rng_uniform_pos(rng1);
        if (r3 <= sigma) {
          Nvec[0]++;
        } else {
          Nvec[r]++;
        }
      }
      N++;
    } else if (r < 6) { // death
      if (r == 3) {
        Nvec[0]--;
      } else if (r == 4) {
        Nvec[1]--;
      } else {
        Nvec[2]--;
      }
      N--;
    } else if (r == 6) { // transfer of a wildtype plasmid
      Nvec[0]--;
      Nvec[1]++;
    } else { // transfer of a mutant plasmid (either from a local cell or an immigrant cell)
      Nvec[0]--;
      Nvec[2]++;
    }

  }
  // return : 0 <-> population extinction, 1 <-> population rescue
  return (N!=0);
}




int main(int argc, char const *argv[]) {

  rescueinput i=filetoinput(argv[1]);
  int NSim = 10000;   // 10000

  double result=0.;
  for (size_t k = 0; k < NSim; k++) {
    bool b=rescue(i.beta, i.phi0, i.gamma, i.sM, i.c, i.sigma, i.alpha, i.u, (int)i.N0, i.s0, k);
    result+= b;
  }
  result/=NSim;

  std::ofstream log("resultsSimulation.txt", std::ios_base::app | std::ios_base::out);
  log << i.beta << " " << i.phi0 << " " << i.gamma << " " << i.sM << " " << i.c << " " << i.sigma << " " << i.alpha << " " << i.u << " " << i.N0 << " " << i.s0 << " " << NSim << " " << result << "\n";

  return 0;
}
